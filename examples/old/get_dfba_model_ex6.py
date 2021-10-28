from cobra.io import read_sbml_model
from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
import multiprocessing
from multiprocessing import Lock


# DfbaModel instance initialized with cobra model
# fba_model = read_sbml_model(
#     join(dirname(__file__), pardir, "sbml-models", "iJR904.xml.gz")
# )


def get_dfba_model(model_dir):
    """Aerobic growth of *S. cerevisiae* on glucose.

    Organism -> Saccharomyces cerevisiae S288C
    Model stored in http://bigg.ucsd.edu/models/iND750
    """
    fba_model = read_sbml_model(model_dir)

    fba_model.solver = "glpk"
    dfba_model = DfbaModel(fba_model)

    # instances of KineticVariable (default initial conditions are 0.0, but
    # can be set here if wanted e.g. Oxygen)
    V = KineticVariable("Volume")
    X = KineticVariable("Biomass")
    Gluc = KineticVariable("Glucose")
    Eth = KineticVariable("Ethanol")
    Glyc = KineticVariable("Glycerol")

    # add kinetic variables to dfba_model
    dfba_model.add_kinetic_variables([V, X, Gluc, Eth, Glyc])

    # instances of ExchangeFlux
    mu = ExchangeFlux("BIOMASS_SC4_bal")
    v_G = ExchangeFlux("EX_glc__D_e")
    v_E = ExchangeFlux("EX_etoh_e")
    v_H = ExchangeFlux("EX_glyc_e")

    # add exchange fluxes to dfba_model
    dfba_model.add_exchange_fluxes([mu, v_G, v_E, v_H])

    # Here add parameters to dfba_model
    params = {"D": 0.044,
              "Vgmax": 8.5}
    Kg = 0.5
    Gin = 100.0
    D = Parameter(list(params.keys())[0], list(params.values())[0])
    Vgmax = Parameter(list(params.keys())[1], list(params.values())[1])

    dfba_model.add_parameters([D, Vgmax])

    # add rhs expressions for kinetic variables in dfba_model
    dfba_model.add_rhs_expression("Volume", D)
    dfba_model.add_rhs_expression("Biomass", mu * X - D * X / V)
    dfba_model.add_rhs_expression("Glucose", v_G * X + D * (Gin - Gluc) / V)
    dfba_model.add_rhs_expression("Ethanol", v_E * X - D * Eth / V)
    dfba_model.add_rhs_expression("Glycerol", v_H * X - D * Glyc / V)

    # add lower/upper bound expressions for exchange fluxes in dfba_model
    # together with expression that must be non-negative for correct evaluation
    # of bounds
    dfba_model.add_exchange_flux_lb("EX_glc__D_e",
                                    Vgmax * (Gluc / (Kg + Gluc)), Gluc)

    # add initial conditions for kinetic variables in dfba_model biomass
    # (gDW/L), metabolites (g/L)
    dfba_model.add_initial_conditions(
        {
            "Volume": 0.5,
            "Biomass": 0.05,
            "Glucose": 10.0,
            "Ethanol": 0.0,
            "Glycerol": 0.0,
        }
    )

    return dfba_model, params


def modifun(dfba_model):
    """Cute little function to modify the dfba model.

    Everything that's needed after setup should go in here.
    """
    # instances of KineticVariable (default initial conditions are 0.0, but
    # can be set here if wanted e.g. Oxygen)
    V = KineticVariable("Volume")
    X = KineticVariable("Biomass")
    Gluc = KineticVariable("Glucose")
    Eth = KineticVariable("Ethanol")
    Glyc = KineticVariable("Glycerol")
    #
    # add kinetic variables to dfba_model
    dfba_model.add_kinetic_variables([V, X, Gluc, Eth, Glyc])

    # instances of ExchangeFlux
    mu = ExchangeFlux("BIOMASS_SC4_bal")
    v_G = ExchangeFlux("EX_glc__D_e")
    v_E = ExchangeFlux("EX_etoh_e")
    v_H = ExchangeFlux("EX_glyc_e")

    # add exchange fluxes to dfba_model
    dfba_model.add_exchange_fluxes([mu, v_G, v_E, v_H])

    # Here add parameters to dfba_model
    params = {"D": 0.044,
              "Vgmax": 8.5}
    Kg = 0.5
    Gin = 100.0
    D = Parameter(list(params.keys())[0], list(params.values())[0])
    Vgmax = Parameter(list(params.keys())[1], list(params.values())[1])

    dfba_model.add_parameters([D, Vgmax])

    # add rhs expressions for kinetic variables in dfba_model
    dfba_model.add_rhs_expression("Volume", D)
    dfba_model.add_rhs_expression("Biomass", mu * X - D * X / V)
    dfba_model.add_rhs_expression("Glucose", v_G * X + D * (Gin - Gluc) / V)
    dfba_model.add_rhs_expression("Ethanol", v_E * X - D * Eth / V)
    dfba_model.add_rhs_expression("Glycerol", v_H * X - D * Glyc / V)

    # add lower/upper bound expressions for exchange fluxes in dfba_model
    # together with expression that must be non-negative for correct evaluation
    # of bounds
    dfba_model.add_exchange_flux_lb("EX_glc__D_e",
                                    Vgmax * (Gluc / (Kg + Gluc)), Gluc)

    # add initial conditions for kinetic variables in dfba_model biomass
    # (gDW/L), metabolites (g/L)
    dfba_model.add_initial_conditions(
        {
            "Volume": 0.5,
            "Biomass": 0.05,
            "Glucose": 10.0,
            "Ethanol": 0.0,
            "Glycerol": 0.0,
        }
    )
    dfba_model.solver_data.set_display("none")


class PicklableDFBAModel:
    """A picklable DFBA model.

    Not pretty, but might work. Needs to override/delegate all methods that are
    expected to be called on the dfba model afterwards.
    """

    def __init__(self, file_, modifunc):
        self.file_ = file_
        self.modifunc = modifunc
        with multiprocessing.Lock():
            self.dfba_model = PicklableDFBAModel.create_dfba_model(file_,
                                                                   modifunc)

    def update_parameters(self, par_dict):      # , *args, **kwargs):
        """Update parameters."""
        self.dfba_model.update_parameters(par_dict)
        # return self.dfba_model()#*args, **kwargs

    def simulate(self, t_start, t_end, t_out, *args, **kwargs):
        """Simulate from the model."""
        concentrations, trajectories = self.dfba_model.simulate(t_start,
                                                                t_end, t_out)
        # return self.dfba_model(*args, **kwargs)
        return concentrations, trajectories

    def __getstate__(self) -> dict:
        """Get a serializable representation of the object for pickling."""
        return {
            'file_': self.file_,
            'modifunc': self.modifunc
            # 'lock': self.lock
        }

    def __setstate__(self, state: dict):
        """Create an object from pickled state."""
        self.file_ = state['file_']
        self.modifunc = state['modifunc']
        with multiprocessing.Lock():
            self.dfba_model = PicklableDFBAModel.create_dfba_model(
                self.file_, self.modifunc)
            # self.dfba_model.solver_data.set_display("none")

    @staticmethod
    def create_dfba_model(file_, modifunc):
        fba_model = read_sbml_model(file_)
        fba_model.solver = "glpk"
        dfba_model = DfbaModel(fba_model)
        modifunc(dfba_model)
        return dfba_model


# import cloudpickle as pickle
# dump = pickle.dumps(obj);
# pickle.loads(dump)
