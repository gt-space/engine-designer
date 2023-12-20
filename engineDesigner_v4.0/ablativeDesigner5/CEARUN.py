from rocketcea.cea_obj import CEA_Obj
ox = "LOX"
fuel = "JetA"
con_rat = 4.514
exp_rat = 2
P_inj_psi= 500
MR= []

ispObj = CEA_Obj(oxName=ox, fuelName=fuel, fac_CR=con_rat)

s = ispObj.get_full_cea_output( Pc=P_inj_psi, MR=MR, eps=exp_rat, short_output=1, pc_units='psia',fac_CR=con_rat)
print(s)