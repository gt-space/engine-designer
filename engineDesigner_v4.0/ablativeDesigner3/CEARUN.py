from rocketcea.cea_obj import CEA_Obj
ox = 'N2O'
fuel = 'Isopropanol'
con_rat = 4.514
exp_rat = 4.5
P_inj_psi= 750
MR= 4

'''thrust = 800*4.44822 # Thrust [N]
P_inj_psi = 18*14.5038 # Chamber pressure at injector face [bar to psi]
exp_rat = 4.5
con_rat = 4.5 #  Contraction ratio sized using correlation in code
MR = 1.8 # Mixture Ratio
fuel = "JetA"
ox="LOX"'''

ispObj = CEA_Obj(oxName=ox, fuelName=fuel, fac_CR=con_rat)

s = ispObj.get_full_cea_output( Pc=P_inj_psi, MR=MR, eps=exp_rat, short_output=1, pc_units='psia',fac_CR=con_rat)
print(s)