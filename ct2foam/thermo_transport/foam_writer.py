import numpy as np
import cantera as ct


def write_species_list(file_name, species_names):
    """
    Writes a species.foam file with an OpenFOAM formatted list of species.
    """
    n_spec=0
    with open(file_name, 'a') as output:
        output.write('species\n(\n')
        for sp_i in species_names:
            n_spec+=1
            output.write('\t')
            output.write(sp_i)
            output.write('\n')
        output.write(');\n\n')

def write_reactions(file_name, input_file):
    """
    Writes an empty reactions.foam file with an OpenFOAM formatted list of 0 reactions.
    """
    gas = ct.Solution(input_file)
    rows = 20
    cols = 20
    n_react=0


    # effi = np.zeros((rows, cols))
    
    
    with open(file_name, 'a') as output:
        # output.write('// - This dictionary intentionally left blank.\n')
        # output.write('reactions\n{\n')
        # output.write('}\n\n')        
        
        output.write(f"reactions\n"+"{\n")

        for r in gas.reactions():
        # Skip falloff / pressure-dependent ones for simplicity
            if str(r.reaction_type)=="reaction":
                n_react += 1
                reactants = " + ".join(f"{v} {k}" for k, v in r.reactants.items())
                products  = " + ".join(f"{v} {k}" for k, v in r.products.items())
                arrhenius = r.rate.pre_exponential_factor, r.rate.temperature_exponent, ((r.rate.activation_energy)/1000)
                output.write(f"    reaction-"+str(n_react)+" \n    {\n")
                output.write(f"        type reversibleArrhenius;\n")
                output.write(f"        reaction \"{reactants} = {products}\";\n")
                output.write(f"        A {arrhenius[0]};\n")
                output.write(f"        beta {arrhenius[1]};\n")
                output.write(f"        Ta {arrhenius[2]};\n")
                output.write("    }\n\n")
        for r in gas.reactions():
            # n_three=0
            effi=[[0 for _ in range(cols)] for _ in range(rows)]
            effi_name=[[0 for _ in range(cols)] for _ in range(rows)]
            effi_=[[0 for _ in range(cols)] for _ in range(rows)]
            if str(r.reaction_type)=="three-body":
                n_react+=1
                reactants = " + ".join(f"{v} {k}" for k, v in r.reactants.items())
                products  = " + ".join(f"{v} {k}" for k, v in r.products.items())
                arrhenius = r.rate.pre_exponential_factor, r.rate.temperature_exponent, ((r.rate.activation_energy)/1000)

                
                # n_three+=1
                # effi=[[0 for _ in range(cols)] for _ in range(rows)]
                # effi_name = [[0 for _ in range(cols)] for _ in range(rows)]
                # effi_=[[0 for _ in range(cols)] for _ in range(rows)]
                
                
                re = str(r.efficiencies)
                x1 = re.replace("{","")
                x2 = x1.replace("}","")

                # effi[(n_three-1)]=x2.split(",")
                effi=x2.split(",")
                
                for i in range(len(effi)):
                    # for j in range(len(effi[i])):
                    ij_elm=str(effi[i])


                    ij_elm=  ((ij_elm.replace(":","")).replace(" '","")).replace("'","")
                    effi_[i]=ij_elm.split(" ")


                spec = [[0 for _ in range(2)] for _ in range(len(gas.species()))]
                n_spec = 0
                for s in gas.species():
                    s=str(s)
                    spec[n_spec][0] = (s.replace("<Species ","")).replace(">","")
                    spec[n_spec][1] = '1'
                    n_spec+=1
                spec.sort()

                effi_name = [row for row in effi_ if not (isinstance(row, list) and all(v == 0 for v in row))]

                print(effi_name)
                print(spec)
                for i in range(len(spec)):
                    for s in range(len(effi_name)):
                        if spec[i][0]==effi_name[s][0]:
                            spec[i][1]=effi_name[s][1]
                        else:
                            continue
                print(spec)




                output.write(f"    reaction-"+str(n_react)+" \n    {\n")
                output.write(f"        type reversibleThirdBodyArrhenius;\n")
                output.write(f"        reaction \"{reactants} = {products}\";\n")
                output.write(f"        A {arrhenius[0]};\n")
                output.write(f"        beta {arrhenius[1]};\n")
                output.write(f"        Ta {arrhenius[2]};\n")

                output.write(f"        coeffs\n")
                output.write(f""+str(len(gas.species()))+"\n(\n")

            
                for i in spec:
                    output.write(f"(" + str(i[0]) + " " +str(i[1]) + ")\n")

                output.write(f")\n"+";\n")

                output.write("    }\n\n")
        output.write(f""+"}")




def write_thermo_transport(file_name, name, MW, As, Ts, poly_mu, poly_kappa, logpoly_mu, logpoly_kappa, nasa7_Tmid, nasa7_Tlo, nasa7_Thi, nasa7_lo, nasa7_hi, elements=None):
    """
    Writes thermophysicalProperties file required dictionary entries according to given data.
    file_name: output file name.
    name: species name / mixture name
    MW: molecular weight
    As, Ts: Sutherland entries
    poly*: transport polynomial entries
    nasa7*: NASA7 polynomial coefficient entry data
    elements: a dictionary with the following syntax:  elements = {"C": 1, "H":1}
    """
    # Numpy-based polynomial fit has a reversed order to OpenFoam dictionary definition
    poly_mu_rev = np.copy(poly_mu)
    poly_mu_rev = np.flip(poly_mu_rev)
    logpoly_mu_rev = np.copy(logpoly_mu)
    logpoly_mu_rev = np.flip(logpoly_mu_rev)

    poly_kappa_rev = np.copy(poly_kappa)
    poly_kappa_rev = np.flip(poly_kappa_rev)
    logpoly_kappa_rev = np.copy(logpoly_kappa)
    logpoly_kappa_rev = np.flip(logpoly_kappa_rev)

    #Write the thermo output in openfoam format
    with open(file_name,'a') as output:
        
        output.write(name+'\n{\n')
        output.write('\t')
        ############################################################################# 
        output.write('specie\n\t{\n')
        output.write('\t\tnMoles \t 1;\n')
        output.write('\t\tmolWeight \t'+repr(MW)+';')
        output.write('\n\t}\n\n')
        #############################################################################

        ############################################################################# 
        output.write('\tthermodynamics')
        #############################################################################         
        output.write('\n\t{\n')
        output.write('\t\tTlow\t\t'+repr(nasa7_Tlo)+';\n')
        output.write('\t\tThigh\t\t'+repr(nasa7_Thi)+';\n')       
        output.write('\t\tTcommon\t\t'+repr(nasa7_Tmid)+';\n')
        output.write('\t\tlowCpCoeffs\t(\t' )
        for wi in range(0,7): #NASA pol has 7 coeffs
            output.write(repr(nasa7_lo[wi])) 
            output.write(' ')
        output.write(' );\n')
        output.write('\t\thighCpCoeffs\t(\t' )
        for wi in range(0,7):#NASA pol has 7 coeffs
            output.write(repr(nasa7_hi[wi])) 
            output.write(' ')
        output.write(' );\n')
        output.write('\t}\n\n')
        
        ############################################################################# 
        output.write('\ttransport \n\t{\n')
        #############################################################################

        output.write('\t\tAs\t'+repr(As)+';\n' )
        output.write('\t\tTs\t'+repr(Ts)+';\n' )

        output.write('\t\tmuLogCoeffs<8>\t(\t' )
        for wi in range(8):
            if(wi < len(logpoly_mu_rev)):
                output.write(repr(logpoly_mu_rev[wi])) 
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')

        output.write('\t\tmuCoeffs<8>\t(\t' )
        for wi in range(8):
            if(wi < len(poly_mu_rev)):
                output.write(repr(poly_mu_rev[wi])) 
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')
        #############################################################################
        output.write('\t\tkappaLogCoeffs<8>\t(\t' )
        for wi in range(8):
            if(wi < len(logpoly_kappa_rev)):
                output.write(repr(logpoly_kappa_rev[wi])) 
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')
        output.write('\t\tkappaCoeffs<8>\t(\t' )
        for wi in range(8):
            if(wi < len(poly_kappa_rev)):
                output.write(repr(poly_kappa_rev[wi])) 
            else:
                output.write("0") 
            output.write(' ')
        output.write(' );\n')
        output.write('\t}\n\n')
        ############################################################################# 

        ############################################################################# 
        if(elements is not None):
            output.write('\telements')
            #############################################################################      
            output.write('\n\t{\n')   
            # Variables for the elemental composition entry:
            for elem_i in elements.keys():
                output.write("\t\t" + elem_i+"\t"+repr(int(elements[elem_i]))+';\n')         
            output.write('\t}\n')
            #############################################################################
        output.write('}\n\n')    

