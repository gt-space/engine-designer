import cantera as ct
# from Cantera_Properties import CanteraProperties

LOX = ct.Solution('Cantera YAML Files/LOX.yaml','liquid_oxygen')
JetA = ct.Solution('Cantera YAML Files/JetA.yaml','jet_a')

# extract all species in the NASA database
full_species = {S.name: S for S in ct.Species.list_from_file('nasa_gas.yaml')}

for spec in full_species:
    print(spec)

# extract only the relevant species
species = [full_species[S] for S in (
    'CH4', 'CO', 'CO2', 'C2H4', 'H2', 'H2O', 'C'
    )]
gas = ct.Solution(thermo='ideal-gas', species=species)

MR_weight = 1.8
MR = MR_weight * 166 / 32

# self.gas = ct.Solution('Cantera YAML Files/CombProducts.yaml')
mixture = ct.Mixture([(LOX,MR/(1+MR)),(JetA,1/(1+MR)),(gas,0)]) # does quantity matter?