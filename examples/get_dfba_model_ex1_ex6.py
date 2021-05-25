from cobra.io import read_sbml_model
from dfba import DfbaModel, ExchangeFlux, KineticVariable, Parameter
import multiprocessing
from multiprocessing import Lock


# DfbaModel instance initialized with cobra model
# fba_model = read_sbml_model(
#     join(dirname(__file__), pardir, "sbml-models", "iJR904.xml.gz")
# )


def setup_dfba_model_ex1(dfba_model):
    """
    Set up the DFBA model with defining rhs, exchange fluxes, parameters,
    initial conditions.
    Example 1.
    Anaerobic growth of *E. coli* on glucose and xylose.

    Organism -> _Escherichia coli str. K-12 substr. MG1655_
    Model in http://bigg.ucsd.edu/models/iJR904
    """
    print("Example 1.")
    # instances of KineticVariable (default initial conditions are 0.0, but
    # can be set here if wanted e.g. Oxygen)
    X = KineticVariable("Biomass")
    Gluc = KineticVariable("Glucose")
    Xyl = KineticVariable("Xylose")
    Oxy = KineticVariable("Oxygen", initial_condition=0.24)
    Eth = KineticVariable("Ethanol")

    # add kinetic variables to dfba_model
    dfba_model.add_kinetic_variables([X, Gluc, Xyl, Oxy, Eth])

    # instances of ExchangeFlux
    mu = ExchangeFlux("BiomassEcoli")
    v_G = ExchangeFlux("EX_glc(e)")
    v_Z = ExchangeFlux("EX_xyl_D(e)")
    v_O = ExchangeFlux("EX_o2(e)")
    v_E = ExchangeFlux("EX_etoh(e)")

    # add exchange fluxes to dfba_model
    dfba_model.add_exchange_fluxes([mu, v_G, v_Z, v_O, v_E])

    # Here add parameters to dfba_model
    params = {"K_g": 0.0027,
              "v_gmax": 10.5,
              "K_z": 0.0165,
              "v_zmax": 6.0}
    K_g = Parameter(list(params.keys())[0], list(params.values())[0])
    v_gmax = Parameter(list(params.keys())[1], list(params.values())[1])
    K_z = Parameter(list(params.keys())[2], list(params.values())[2])
    v_zmax = Parameter(list(params.keys())[3], list(params.values())[3])

    dfba_model.add_parameters([K_g, v_gmax, K_z, v_zmax])
    # dfba_model.add_parameters([K_g,K_z])#,init_biomass,init_glucose,
                              # init_xylose,init_oxygen,init_ethanol])

    # add rhs expressions for kinetic variables in dfba_model
    dfba_model.add_rhs_expression("Biomass", mu * X)
    dfba_model.add_rhs_expression("Glucose", v_G * 180.1559 * X / 1000.0)
    dfba_model.add_rhs_expression("Xylose", v_Z * 150.13 * X / 1000.0)
    dfba_model.add_rhs_expression("Oxygen", v_O * 16.0 * X / 1000.0)
    dfba_model.add_rhs_expression("Ethanol", v_E * 46.06844 * X / 1000.0)

    # add lower/upper bound expressions for exchange fluxes in dfba_model
    # together with expression that must be non-negative for correct evaluation
    # of bounds
    # v_gmax = 10.5
    # K_g = 0.0027
    dfba_model.add_exchange_flux_lb(
        "EX_glc(e)", v_gmax * (Gluc / (K_g + Gluc)) * (1 / (1 + Eth / 20.0)),
        Gluc)   # v_g glucose
    dfba_model.add_exchange_flux_lb("EX_o2(e)", 15.0 * (Oxy / (0.024 + Oxy)),
                                    Oxy)    # v_o, oxygen
    # v_zmax = 6.0
    # K_z = 0.0165
    dfba_model.add_exchange_flux_lb(
        "EX_xyl_D(e)",
        v_zmax * (Xyl / (K_z + Xyl)) * (1 / (1 + Eth / 20.0)) *
        (1 / (1 + Gluc / 0.005)), Xyl)  # v_z, xylose

    # add initial conditions for kinetic variables in dfba_model biomass
    # (gDW/L), metabolites (g/L)
    dfba_model.add_initial_conditions(
        {
            "Biomass": 0.03,
            "Glucose": 15.5,
            "Xylose": 8.0,
            "Oxygen": 0.0,
            "Ethanol": 0.0,
        }
    )
    return dfba_model, params


def setup_dfba_model_ex1_aerobic(dfba_model):
    """
    Set up the DFBA model with defining rhs, exchange fluxes, parameters,
    initial conditions.
    Example 1 - Real Data, Eitemann et al 2008.
    Aerobic growth of *E. coli* on glucose and xylose.

    Organism -> _Escherichia coli str. K-12 substr. MG1655_
    Model in http://bigg.ucsd.edu/models/iJR904
    """
    print("Example 1. Real Data. Aerobic Growth.")
    # instances of KineticVariable (default initial conditions are 0.0, but
    # can be set here if wanted e.g. Oxygen)
    X = KineticVariable("Biomass")
    Gluc = KineticVariable("Glucose")
    Xyl = KineticVariable("Xylose")
    Eth = KineticVariable("Ethanol")

    # add kinetic variables to dfba_model
    dfba_model.add_kinetic_variables([X, Gluc, Xyl, Eth])

    # instances of ExchangeFlux
    mu = ExchangeFlux("BiomassEcoli")
    v_G = ExchangeFlux("EX_glc(e)")
    v_Z = ExchangeFlux("EX_xyl_D(e)")
    v_E = ExchangeFlux("EX_etoh(e)")

    # add exchange fluxes to dfba_model
    dfba_model.add_exchange_fluxes([mu, v_G, v_Z, v_E])

    # Here add parameters to dfba_model
    params = {"K_g": 0.0027,
              "v_gmax": 10.5,
              "K_z": 0.0165,
              "v_zmax": 6.0}
    K_g = Parameter(list(params.keys())[0], list(params.values())[0])
    v_gmax = Parameter(list(params.keys())[1], list(params.values())[1])
    K_z = Parameter(list(params.keys())[2], list(params.values())[2])
    v_zmax = Parameter(list(params.keys())[3], list(params.values())[3])

    dfba_model.add_parameters([K_g, v_gmax, K_z, v_zmax])

    # add rhs expressions for kinetic variables in dfba_model
    dfba_model.add_rhs_expression("Biomass", mu * X)
    dfba_model.add_rhs_expression("Glucose", v_G * 180.1559 * X / 1000.0)
    dfba_model.add_rhs_expression("Xylose", v_Z * 150.13 * X / 1000.0)
    dfba_model.add_rhs_expression("Ethanol", v_E * 46.06844 * X / 1000.0)

    # add lower/upper bound expressions for exchange fluxes in dfba_model
    # together with expression that must be non-negative for correct evaluation
    # of bounds
    # v_gmax = 10.5
    # K_g = 0.0027
    dfba_model.add_exchange_flux_lb(
        "EX_glc(e)", v_gmax * (Gluc / (K_g + Gluc)) * (1 / (1 + Eth / 20.0)),
        Gluc)  # v_g glucose
    # v_zmax = 6.0
    # K_z = 0.0165
    dfba_model.add_exchange_flux_lb(
        "EX_xyl_D(e)",
        v_zmax * (Xyl / (K_z + Xyl)) * (1 / (1 + Eth / 20.0)) *
        (1 / (1 + Gluc / 0.005)), Xyl)  # v_z, xylose

    # add initial conditions for kinetic variables in dfba_model biomass
    # (gDW/L), metabolites (g/L)
    dfba_model.add_initial_conditions(
        {
            "Biomass": 0.03,
            "Glucose": 15.5,
            "Xylose": 8.0,
            "Ethanol": 0.0,
        }
    )
    return dfba_model, params


def setup_dfba_model_ex6(dfba_model):
    """
    Set up the DFBA model with defining rhs, exchange fluxes, parameters,
    initial conditions.
    Example 6.
    Aerobic growth of *S. cerevisiae* on glucose.

    Organism -> Saccharomyces cerevisiae S288C
    Model stored in http://bigg.ucsd.edu/models/iND750
    """
    print("Example 6.")
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


def get_dfba_model(model_dir, example_name):
    """
    Read FBA model from SBML and build DFBA model with defining parameters.

    Parameters
    ----------
    model_dir: str
        model directory to sbml-FBA-model
    example_name: str
        name of dfba-example (example1 | example6)
    Returns
    -------
    dfba_model
    params: dict
        Dictionary with Parameter names as Keys, and their parameter values
    """
    # TODO: this function can be optimized/should be rewritten, if its only
    #  used for getting param_dict, then the model does not need to be loaded.
    fba_model = read_sbml_model(model_dir)
    if example_name=="example1_aerobic":
        fba_model.reactions.get_by_id("EX_o2(e)").lower_bound = - 15.0 * (
                    0.24 / (0.024 + 0.24))

    fba_model.solver = "glpk"
    dfba_model = DfbaModel(fba_model)

    if example_name=="example1":
        dfba_model, params = setup_dfba_model_ex1(dfba_model)
    elif example_name=="example6":
        dfba_model, params = setup_dfba_model_ex6(dfba_model)
    elif example_name=="example1_aerobic":
        dfba_model, params = setup_dfba_model_ex1_aerobic(dfba_model)
    else:
        raise ValueError("Choose either 'example1' or 'example1_aerobic' "
                         "or 'example6' as argument for examplename.")

    return dfba_model, params


def modifun(dfba_model, example_name):
    """Cute little function to modify the dfba model.

    Everything that's needed after setup should go in here.
    """
    if example_name == "example1":
        setup_dfba_model_ex1(dfba_model)
    elif example_name == "example6":
        setup_dfba_model_ex6(dfba_model)
    elif example_name == "example1_aerobic":
        dfba_model, params = setup_dfba_model_ex1_aerobic(dfba_model)
    else:
        raise ValueError("Choose either 'example1' or 'example1_aerobic' "
                         "or 'example6' as argument for examplename.")

    dfba_model.solver_data.set_display("none")


class PicklableDFBAModel:
    """A picklable DFBA model.

    Not pretty, but might work. Needs to override/delegate all methods that are
    expected to be called on the dfba model afterwards.
    """

    def __init__(self, file_, modifunc, example_name):
        self.file_ = file_
        self.modifunc = modifunc
        self.example_name = example_name
        with multiprocessing.Lock():
            self.dfba_model = PicklableDFBAModel.create_dfba_model(file_,
                                                                   modifunc,
                                                                   example_name)

    def update_parameters(self, par_dict):#, *args, **kwargs):
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
            'modifunc': self.modifunc,
            'example_name': self.example_name
            # 'lock': self.lock
        }

    def __setstate__(self, state: dict):
        """Create an object from pickled state."""
        self.file_ = state['file_']
        self.modifunc = state['modifunc']
        self.example_name = state['example_name']
        with multiprocessing.Lock():
            self.dfba_model = PicklableDFBAModel.create_dfba_model(
                self.file_, self.modifunc, self.example_name)

    @staticmethod
    def create_dfba_model(file_, modifunc, example_name):
        fba_model = read_sbml_model(file_)
        if example_name == "example1_aerobic":
            fba_model.reactions.get_by_id("EX_o2(e)").lower_bound = - 15.0 * (
                    0.24 / (0.024 + 0.24))
        fba_model.solver = "glpk"
        dfba_model = DfbaModel(fba_model)
        modifunc(dfba_model, example_name)
        return dfba_model


# import cloudpickle as pickle
# dump = pickle.dumps(obj);
# pickle.loads(dump)
