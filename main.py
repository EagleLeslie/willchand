import os
import numpy as np

def new_intf_file(intc):
	"""
	Loop through new interface files
	"""
	os.system(f'rm INTFC')
	os.system(f'touch INTFC')
	with open('INTFC','r') as file:
		lines = file.readlines()

	lines = f"{intc}"

	with open('INTFC','w') as file:
		file.writelines(lines)

if __name__ == '__main__':
	
	# atoms in system:
	fe = 68
	mg = 16
	o = 48
	si = 16

	vol = 1356.7902996719795

	interface_names = ["Fe", "Mg", "O", "Si"]
	interface_consts = np.array([fe/vol, mg/vol, o/vol, si/vol])

	for i,x in enumerate(interface_consts):
		intc = x
		new_intf_file(intc)
		print("************************************************************\n")
		print(f"INTERFACE CONSTANT OF ATOM {interface_names[i]} IS {intc} \n")
		print("************************************************************\n")

		os.system("./main")
		os.system(f"mv intf1 intf1-{interface_names[i]}")
		os.system(f"mv intf2 intf2-{interface_names[i]}")

