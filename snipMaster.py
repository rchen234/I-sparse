import numpy as np
import ast
from gurobipy import *
import master
import time

class SNIPinst(master.StochIPinst):
	def __init__(self,instance,R,snipNo):
		self.instance = instance
		self.name = str(instance)+'_'+str(R)+'_'+str(snipNo)
		self.R = R
		self.snipNo = snipNo
		self.Nscen = None
		self.Nfv = None
		self.cVec = None
		self.pVec = None
		self.BendersCuts = {}
		self.A = []
		self.AD = []
		self.Alist = None
		self.ADlist = None
		self.ALL = []
		self.N = []
		self.SCEN = []
		self.PSI = []
		self.ind = {}
		self.n = None
		self.na = None
		self.a = None
		self.ad = None
		self.PrimalMaster = Model('Master')
		self.theta = {}
		self.x = {}
		# Benders cuts lists
		self.cutlist = {}
		self.NonzeroInd = {}
		self.NonzeroSparse = {}
		# All cuts coef. lists
		self.coeflist = {}
		self.thetaCutList = {}
		self.Ncuts = 0
		self.Nsubs = 0
		# Subproblems
		self.scensub = {}
		self.scenrc = None
		self.subprob = {}
		# Stored solutions
		self.Solns = {}
		# Graph edges
		self.Edges = {}
		# Time
		self.GeneralSparseT = 0.0
		self.BendersSparseT = 0.0
		self.MixedSparseT = 0.0
		self.RhsT = 0.0
		self.ExactlyEvaluated = {}
		self.IEvaluated = {}
		self.Iineqs = {}
		self.TiltedIneqs = {}

	def Initialize(self):
		objCoef = {}
		for i in range(self.Nfv):
			objCoef[i] = self.cVec[i]
		regCoefy = 1.0
		t0 = time.time()
		for s in range(self.Nscen):
			self.scensub[s] = Model()
			self.scensub[s].setParam('OutputFlag', False)
			self.scensub[s].modelSense = GRB.MINIMIZE
			x_scen = {}
			for i in range(self.Nfv):
				x_scen[i] = self.scensub[s].addVar(vtype=GRB.BINARY,obj=objCoef[i],name="x"+str(i))
			ypi = {}
			for i in range(self.n):
				ypi[i] = self.scensub[s].addVar(lb=0.0,ub=1.0,obj=0.0,name="y"+str(i))
			ypi[self.ind[self.SCEN[s,0]]].obj = regCoefy
			self.scensub[s].addConstr(quicksum(x_scen[i] for i in range(self.Nfv)) <= self.R)
			for i in range(self.a):
				self.scensub[s].addConstr(ypi[self.ind[self.A[i][0]]] - self.A[i][2]*ypi[self.ind[self.A[i][1]]] >= 0)
			for i in range(self.ad):
				self.scensub[s].addConstr(ypi[self.ind[self.AD[i][0]]] - self.AD[i][2]*ypi[self.ind[self.AD[i][1]]] + (self.AD[i][2] - self.AD[i][3])*self.PSI[s,self.ind[self.AD[i][1]]]*x_scen[i] >= 0,name="con"+str(i))
				self.scensub[s].addConstr(ypi[self.ind[self.AD[i][0]]] - self.AD[i][3]*ypi[self.ind[self.AD[i][1]]] >= 0)
			self.scensub[s].addConstr(ypi[self.ind[self.SCEN[s,1]]] == 1,name="cont")
			self.scensub[s].setParam( 'OutputFlag', False )
			self.scensub[s].modelSense = GRB.MINIMIZE
			self.scensub[s].update()
			
		print("Build scenario subproblems: "+str(time.time()-t0))
		t0 = time.time()
		self.scenrc = Model()
		self.scenrc.setParam('OutputFlag', False)
		self.scenrc.modelSense = GRB.MINIMIZE
		ypi = {}
		subgConstr = {}
		for i in range(self.n):
			ypi[i] = self.scenrc.addVar(lb=0.0,ub=1.0,obj=0.0,name="y"+str(i))
		for i in range(self.a):
			self.scenrc.addConstr(ypi[self.ind[self.A[i][0]]] - self.A[i][2]*ypi[self.ind[self.A[i][1]]] >= 0)
		for i in range(self.ad):
			subgConstr[i] = self.scenrc.addConstr(ypi[self.ind[self.AD[i][0]]] - self.AD[i][2]*ypi[self.ind[self.AD[i][1]]] >= 0,name="subgCon"+str(i))
			self.scenrc.addConstr(ypi[self.ind[self.AD[i][0]]] - self.AD[i][3]*ypi[self.ind[self.AD[i][1]]] >= 0)
		self.scenrc.update()
		print("Build recourse subproblems: "+str(time.time()-t0))

	def FreeMemory(self):
		for s in range(self.Nscen):
			del self.scensub[s]
		del self.scenrc
		del self.PrimalMaster

	def SolveExtensive(self):
		Extensive = Model('Extensive')
		Extensive.modelSense = GRB.MINIMIZE
		x = {}
		ypi = {}
		for i in range(self.Nfv):
			x[i] = Extensive.addVar(vtype=GRB.BINARY)
		for s in range(self.Nscen):
			for i in range(self.n):
				ypi[s,i] = Extensive.addVar(lb=0.0)
		Extensive.addConstr(quicksum(x[i] for i in range(self.Nfv)) <= self.R)
		for s in range(self.Nscen):
			for i in range(self.a):
				Extensive.addConstr(ypi[s,self.ind[self.A[i][0]]] - self.A[i][2]*ypi[s,self.ind[self.A[i][1]]] >= 0)
			for i in range(self.ad):
				Extensive.addConstr(ypi[s,self.ind[self.AD[i][0]]] - self.AD[i][2]*ypi[s,self.ind[self.AD[i][1]]] >= -(self.AD[i][2] - self.AD[i][3])*self.PSI[s,self.ind[self.AD[i][1]]]*x[i])
				Extensive.addConstr(ypi[s,self.ind[self.AD[i][0]]] - self.AD[i][3]*ypi[s,self.ind[self.AD[i][1]]] >= 0)
			Extensive.addConstr(ypi[s,self.ind[self.SCEN[s,1]]] == 1)
		Extensive.setObjective(quicksum(self.pVec[s]*ypi[s,self.ind[self.SCEN[s,0]]] for s in range(self.Nscen)))
		Extensive.update()
		t0 = time.time()
		Extensive.optimize()
		print('Solve Extensive Form: '+str(self.name))
		print('Optimal Obj. Value: '+str(Extensive.objval)+', Total Time: '+str(time.time()-t0))



	def BuildPrimal(self):
		for s in range(self.Nscen):
			self.theta[s] = self.PrimalMaster.addVar(lb=-GRB.INFINITY,name="theta"+str(s))
		for i in range(self.Nfv):
			self.x[i] = self.PrimalMaster.addVar(lb=0.0,ub=1.0,name="x"+str(i))
		self.PrimalMaster.addConstr(quicksum(self.x[i] for i in range(self.Nfv)) <= self.R)
		self.PrimalMaster.setObjective(quicksum(self.pVec[s]*self.theta[s] for s in range(self.Nscen))+quicksum(self.cVec[j]*self.x[j] for j in range(self.Nfv)))
		self.PrimalMaster.setParam( 'OutputFlag', False )
		self.PrimalMaster.modelSense = GRB.MINIMIZE
		self.PrimalMaster.update()
		self.addBenders(init=True)

	def readData(self):
		instanceName = self.instance+'.txt'
		arcFile = "Instances/SNIP/nonint_"+instanceName
		intArcFile = "Instances/SNIP/int_"+instanceName
		f = open(arcFile)
		x = f.readlines()
		f.close()
		for line in x:
			if x.index(line)%2 == 0:
				z = line.strip('\n').split('\t\t')
				z = [int(z[0]),int(z[1]),float(z[2]),float(z[2])]
				self.N.append(int(z[0]))
				self.N.append(int(z[1]))
				# Unsorted
				self.A.append(z)

		f = open(intArcFile)
		x = f.readlines()
		f.close()
		for line in x:
			if x.index(line)%2 == 0:
				z = line.strip('\n').split('\t\t')
				if self.snipNo == 1:
					z = [int(z[0]),int(z[1]),float(z[2]),float(z[3])]
				elif self.snipNo == 2:
					z = [int(z[0]),int(z[1]),float(z[2]),0.5*float(z[2])]
				elif self.snipNo == 3:
					z = [int(z[0]),int(z[1]),float(z[2]),0.1*float(z[2])]
				elif self.snipNo == 4:
					z = [int(z[0]),int(z[1]),float(z[2]),0.0]
				# Unsorted
				self.N.append(int(z[0]))
				self.N.append(int(z[1]))
				self.AD.append(z)

		self.N = np.sort(np.unique(np.array(self.N)))
		self.ALL = self.A+self.AD

		self.a = np.size(self.A,0)
		self.ad = np.size(self.AD,0)
		self.Nfv = self.ad
		self.cVec = [0.0]*self.Nfv
		self.na = np.size(self.ALL,0)
		self.n = np.size(self.N,0)

		f = open("Instances/SNIP/Scenarios.txt")
		x = f.readlines()
		f.close()
		for line in x:
			z = line.strip('\n').split('\t')
			z = [int(z[0]),int(z[1]),float(z[2])]
			self.SCEN.append(z)
		self.SCEN = np.array(self.SCEN)
		self.pVec = self.SCEN[:,2]

		f = open("Instances/SNIP/psi_reformat.txt")
		x = f.readlines()
		f.close()

		for string in x:
			z = string.strip('\n').split(' ')
			z = [float(i) for i in z]
			self.PSI.append(z)
		self.PSI = np.array(self.PSI)

		for i in range(self.n):
			self.ind[self.N[i]] = i

		ns = np.size(self.SCEN,0)
		self.Nscen = ns
		for s in range(self.Nscen):
			self.cutlist[s] = []
			self.coeflist[s] = []
			self.thetaCutList[s] = []
			self.BendersCuts[s] = np.empty((0,self.Nfv+1), float)
			self.NonzeroInd[s] = []
			self.NonzeroSparse[s] = []
			self.Solns[s] = []
			self.ExactlyEvaluated[s] = {}
			self.IEvaluated[s] = []
			self.Iineqs[s] = []
			self.TiltedIneqs[s] = []

		self.Alist = [(self.ind[self.A[i][0]],self.ind[self.A[i][1]]) for i in range(self.a)]
		self.ADlist = [(self.ind[self.AD[i][0]],self.ind[self.AD[i][1]]) for i in range(self.ad)]

		self.Initialize()


	
	def addGeneralSparse(self,s,x_input,I,shift=[0]):
		t0 = time.time()
		I.sort()
		IsubList = []
		for card in range(min(len(I)+1,self.R)):
			IsubList = IsubList+list(itertools.combinations(I,card))
		SepProb = Model()
		mu = SepProb.addVars(I, lb=-GRB.INFINITY, name="mu")
		gamma = SepProb.addVar(lb=-GRB.INFINITY, name="gamma")
		SepProb.setObjective(quicksum(x_input[i]*mu[i] for i in I)+gamma)
		SepProb.modelSense = GRB.MAXIMIZE
		SepProb.setParam( 'OutputFlag', False )
		Ibar = [i for i in range(self.Nfv) if i not in I]
		solns = {}
		for subI in IsubList:
			x_dict = {}
			for i in I:
				if i in subI:
					x_dict[i] = 1
				else:
					x_dict[i] = 0
			if s != -1:
				if (tuple(I),subI) in self.ExactlyEvaluated[s].keys():
					objx = self.ExactlyEvaluated[s][tuple(I),subI]
				else:
					subprob = self.scensub[s].relax()
					subprob.setObjective(subprob.getVarByName("y"+str(self.ind[self.SCEN[s,0]])))
					for i in I:
						if i in subI:
							subprob.getVarByName("x"+str(i)).lb = 1
						else:
							subprob.getVarByName("x"+str(i)).ub = 0
					subprob.update()
					rhsT0 = time.time()
					subprob.optimize()
					self.RhsT += time.time()-rhsT0
					objx = subprob.objval
					self.ExactlyEvaluated[s][tuple(I),subI] = objx
			else:
				objx = 0.0
				for scen in range(self.Nscen):
					if (tuple(I),subI) in self.ExactlyEvaluated[scen].keys():
						objx += self.ExactlyEvaluated[scen][tuple(I),subI]*self.pVec[scen]
					else:
						subprob = self.scensub[scen].relax()
						subprob.setObjective(subprob.getVarByName("y"+str(self.ind[self.SCEN[scen,0]])))
						for i in I:
							if i in subI:
								subprob.getVarByName("x"+str(i)).lb = 1
							else:
								subprob.getVarByName("x"+str(i)).ub = 0
						subprob.update()
						rhsT0 = time.time()
						subprob.optimize()
						self.RhsT += time.time()-rhsT0
						obj = subprob.objval
						self.ExactlyEvaluated[scen][tuple(I),subI] = obj
						objx += obj*self.pVec[scen]
			SepProb.addConstr(quicksum(x_dict[i]*mu[i] for i in I)+gamma<=objx)
		SepProb.update()
		SepProb.optimize()
		mu_soln = {}
		for i in I:
			mu_soln[i] = mu[i].x
		gamma_soln = gamma.x
		self.GeneralSparseT += time.time()-t0
		return SepProb.objval,mu_soln,gamma_soln,solns

	def SolveBendersSub(self,scen_id,x_input):
		for i in range(self.ad):
			self.scenrc.getConstrByName("subgCon"+str(i)).RHS = -(self.AD[i][2] - self.AD[i][3])*self.PSI[scen_id,self.ind[self.AD[i][1]]]*x_input[i]
		self.scenrc.getVarByName("y"+str(self.ind[self.SCEN[scen_id,1]])).lb = 1.0
		self.scenrc.getVarByName("y"+str(self.ind[self.SCEN[scen_id,0]])).obj = 1.0
		self.scenrc.update()
		self.scenrc.optimize()
		ObjValue = self.scenrc.objval
		subg = {}
		#Benderscoef = []
		for i in range(self.ad):
			subg[i] = -(self.AD[i][2] - self.AD[i][3])*self.PSI[scen_id,self.ind[self.AD[i][1]]]*self.scenrc.getConstrByName("subgCon"+str(i)).pi
			#Benderscoef.append(-(self.AD[i][2] - self.AD[i][3])*self.PSI[scen_id,self.ind[self.AD[i][1]]]*self.scenrc.getConstrByName("subgCon"+str(i)).pi)
		const = self.scenrc.objval - sum(subg[i]*x_input[i] for i in range(self.Nfv))
		self.scenrc.getVarByName("y"+str(self.ind[self.SCEN[scen_id,1]])).lb = 0.0
		self.scenrc.getVarByName("y"+str(self.ind[self.SCEN[scen_id,0]])).obj = 0.0
		self.scenrc.update()
		#Benderscoef.append(const)
		Benderscoef = []
		for i in range(self.Nfv):
			Benderscoef.append(subg[i])
		Benderscoef.append(const)
		Benderscoef = np.array([Benderscoef])
		self.BendersCuts[scen_id] = np.append(self.BendersCuts[scen_id], Benderscoef, axis=0)
		return ObjValue,const,subg
	
	def SolveExtensive(self,tLimit=60*60):
		t0 = time.time()
		Ext = Model()
		Ext.Params.LogFile="Results/"+str(self.name)+"_Ext.log"
		x = Ext.addVars(self.Nfv, vtype=GRB.BINARY)
		ypi = Ext.addVars(self.Nscen,self.n)
		Ext.addConstr(quicksum(x[i] for i in range(self.Nfv)) <= self.R)
		Ext.addConstrs(ypi[s,self.ind[self.A[i][0]]] - self.A[i][2]*ypi[s,self.ind[self.A[i][1]]] >= 0 for s in range(self.Nscen) for i in range(self.a))
		Ext.addConstrs(ypi[s,self.ind[self.AD[i][0]]] - self.AD[i][2]*ypi[s,self.ind[self.AD[i][1]]] + (self.AD[i][2] - self.AD[i][3])*self.PSI[s,self.ind[self.AD[i][1]]]*x[i] >= 0 for s in range(self.Nscen) for i in range(self.ad))
		Ext.addConstrs(ypi[s,self.ind[self.AD[i][0]]] - self.AD[i][3]*ypi[s,self.ind[self.AD[i][1]]] >= 0 for s in range(self.Nscen) for i in range(self.ad))
		Ext.addConstrs(ypi[s,self.ind[self.SCEN[s,1]]] == 1 for s in range(self.Nscen))
		Ext.setObjective(quicksum(self.pVec[s]*ypi[s,self.ind[self.SCEN[s,0]]] for s in range(self.Nscen)))
		Ext.Params.timeLimit = tLimit
		Ext.update()
		Ext.optimize()

		wrtStr = str(Ext.objval)+"\t"+str(Ext.ObjBound)+"\t"+str(time.time()-t0)+"\t"+str(Ext.NodeCount)
		f = open("Results/"+str(self.name)+"_Ext.txt","a")
		f.write(wrtStr)
		f.close
