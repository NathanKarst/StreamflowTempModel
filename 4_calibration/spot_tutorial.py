


def plot_best_result(dbname, evaluation):
	import pandas as pd
	from matplotlib import pyplot as plt
	import seaborn
	import numpy as np
	df = pd.read_csv(dbname)
	cols = df.columns.tolist()
	sim_start = cols.index('simulation0')
	df_mat = df.as_matrix()
	max_ind = np.nanargmax(df_mat[:,0])
	print df_mat[max_ind,:sim_start]
	best_flow = df_mat[max_ind,sim_start:]
	fig = plt.figure()
	plt.plot(range(len(best_flow)), best_flow)
	plt.plot(range(len(evaluation)), np.array(evaluation))
	fig.savefig('best_simulation_plot.png')


def main():
	import spotpy
	import os
	from spot_setup import spot_setup
	import pickle
	parent_dir = os.path.dirname(os.path.dirname(os.getcwd()))
	from matplotlib import pyplot as plt
	import numpy as np
	import seaborn as sns
	import os
	import pickle
	from datetime import date
	import pandas as pd
	import numpy as np
	#Create samplers for every algorithm:
	group_params = pickle.load( open( os.path.join(parent_dir,'model_data','group_params.p'), "rb" ))
	param_ranges = pickle.load( open( os.path.join(parent_dir,'model_data','param_ranges.p'), "rb" ))
	results=[]
	spot_setup=spot_setup(group_params, param_ranges)
	rep=2000
	dbname = 'ElderSC.csv'
	sampler=spotpy.algorithms.sceua(spot_setup, dbname='ElderSC', dbformat='csv')
	sampler.sample(rep)
	results = sampler.getdata()                       # Load the results
	#spotpy.analyser.plot_regression(results, spot_setup.evaluation()) 
	spotpy.analyser.plot_parametertrace(results)
	evaluation = np.array(spot_setup.evaluation())
	plot_best_result(dbname, evaluation)



if __name__ == '__main__': 
    main()

# sampler=spotpy.algorithms.lhs(spot_setup,   dbname='RosenLHS',   dbformat='csv')
# results.append(sampler.sample(rep))

# sampler=spotpy.algorithms.mle(spot_setup,   dbname='RosenMLE',   dbformat='csv')
# results.append(sampler.sample(rep))

# sampler=spotpy.algorithms.mcmc(spot_setup,  dbname='ElderMCMC',  dbformat='csv')
# results.append(sampler.sample(rep))

# sampler=spotpy.algorithms.sceua(spot_setup, dbname='RosenSCEUA', dbformat='csv')
# results.append(sampler.sample(rep,ngs=4))

# sampler=spotpy.algorithms.sa(spot_setup,    dbname='RosenSA',    dbformat='csv')
# results.append(sampler.sample(rep))

# sampler=spotpy.algorithms.demcz(spot_setup, dbname='RosenDEMCz', dbformat='csv')
# results.append(sampler.sample(rep,nChains=4))

# sampler=spotpy.algorithms.rope(spot_setup,  dbname='RosenROPE',  dbformat='csv')
# results.append(sampler.sample(rep))

def plot_best_result(dbname, evaluation):
	df = pd.read_csv(dbname)
	cols = df.columns
	sim_start = cols.index('simulation0')
	df_mat = df.as_matrix
	max_ind = np.argmax(df_mat[:,0])
	best_flow = df_mat[max_ind,sim_start:]
	fig = plt.figure()
	plt.plot(range(len(best_flow)), best_flow)
	plt.plot(range(len(evaluation)), np.array(evaluation))
	fig.savefig('best_simulation_plot.png')
