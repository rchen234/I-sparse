import numpy as np
import ast
from gurobipy import *
import master
import time
import pickle

class LLAinst(master.StochIPinst):
	def __init__(self,n,kappa,Z,oneOverc,Cksize):
		self.Nscen = kappa
		self.Nfv = n
		self.c = 1/oneOverc
		self.instance = 'lla_'+str(n)+'_'+str(kappa)+'_'+str(Z)+'_'+str(Cksize)
		self.name = self.instance+'_'+str(oneOverc)
		self.cVec = [0]*n
		self.pVec = [1]*kappa
		self.lam = None
		self.w = None
		self.Ck = None
		self.v = None
		self.v0 = None
		self.BendersCuts = {}
		self.PrimalMaster = None
		self.theta = {}
		self.x = {}
		self.scenrc = {}
		self.scensub = {}
		# Benders cuts lists
		self.cutlist = {}
		self.NonzeroInd = {}
		self.NonzeroSparse = {}
		# All cuts coef. lists
		self.coeflist = {}
		self.thetaCutList = {}
		self.Ncuts = 0
		self.Nsubs = 0
		# Stored solutions
		self.Solns = {}
		# Time
		self.GeneralSparseT = 0.0
		self.BendersSparseT = 0.0
		self.MixedSparseT = 0.0
		self.RhsT = 0.0
		self.ExactlyEvaluated = {}
		self.IEvaluated = {}
		self.Iineqs = {}
		self.TiltedIneqs = {}

	def readData(self,singleScen=True):
		filename = 'Instances/LLA/'+self.instance+'.pkl'
		file = open(filename, 'rb')
		[self.lam,self.w,self.Ck,self.v,self.v0] = pickle.load(file)
		file.close()
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
		self.Initialize()

	def Initialize(self):
		t0 = time.time()
		for s in range(self.Nscen):
			self.scenrc[s] = Model()
			self.scenrc[s].setParam('OutputFlag', False)
			self.scenrc[s].modelSense = GRB.MINIMIZE
			zk = {}
			for i in self.Ck[s]:
				zk[i] = self.scenrc[s].addVar(lb=0.0,obj=-self.w[i]*self.lam[s]*self.v[i,s])
			yk = self.scenrc[s].addVar(lb=0.0,obj=0.0)
			self.scenrc[s].addConstr(self.v0[s]*yk+quicksum(self.v[i,s]*zk[i] for i in self.Ck[s]) == 1)
			for i in self.Ck[s]:
				self.scenrc[s].addConstr(self.v0[s]*yk-self.v0[s]*zk[i] <= 0.0, name="con:1-xi"+str(i))
				self.scenrc[s].addConstr(zk[i] <= yk)
				self.scenrc[s].addConstr((self.v0[s]+self.v[i,s])*zk[i] <= 0.0, name="con:xi"+str(i))
			#self.scenrc[s].setParam('DualReductions', False)
			self.scenrc[s].update()
		print("Build recourse subproblems: "+str(time.time()-t0))

		for s in range(self.Nscen):
			self.scensub[s] = Model()
			self.scensub[s].setParam('OutputFlag', False)
			self.scensub[s].modelSense = GRB.MINIMIZE
			x_scen = {}
			zk = {}
			for i in range(self.Nfv):
				x_scen[i] = self.scensub[s].addVar(lb=0.0,ub=1.0,obj=0.0,name="x"+str(i))
			for i in self.Ck[s]:
				zk[i] = self.scensub[s].addVar(lb=0.0,obj=-self.w[i]*self.lam[s]*self.v[i,s])
			yk = self.scensub[s].addVar(lb=0.0,obj=0.0)
			self.scensub[s].addConstr(quicksum(x_scen[i] for i in range(self.Nfv)) <= self.c*self.Nfv)
			self.scensub[s].addConstr(self.v0[s]*yk+quicksum(self.v[i,s]*zk[i] for i in self.Ck[s]) == 1)
			for i in self.Ck[s]:
				self.scensub[s].addConstr(self.v0[s]*yk-self.v0[s]*zk[i] <= 1-x_scen[i])
				self.scensub[s].addConstr(zk[i] <= yk)
				self.scensub[s].addConstr((self.v0[s]+self.v[i,s])*zk[i] <= x_scen[i])
			self.scensub[s].update()
		print("Build scenario subproblems: "+str(time.time()-t0))

	def FreeMemory(self):
		for s in range(self.Nscen):
			del self.scensub[s]
			del self.scenrc[s]
		del self.PrimalMaster

	def SolveExtensive(self,tLimit=60*60):
		t0 = time.time()
		Ext = Model()
		Ext.Params.LogFile="Results/"+str(self.name)+"_Ext.log"
		x = Ext.addVars(range(self.Nfv), vtype=GRB.BINARY)
		y = Ext.addVars(range(self.Nscen), lb=0.0)
		z = Ext.addVars(self.Nfv, self.Nscen, lb=0.0)
		Ext.setObjective(quicksum(-self.w[i]*self.lam[k]*self.v[i,k]*z[i,k] for k in range(self.Nscen) for i in self.Ck[k]))
		Ext.addConstr(quicksum(x[i] for i in range(self.Nfv)) <= self.c*self.Nfv)
		for k in range(self.Nscen):
			Ext.addConstr(self.v0[k]*y[k]+quicksum(self.v[i,k]*z[i,k] for i in self.Ck[k]) == 1)
			for i in self.Ck[k]:
				Ext.addConstr(self.v0[k]*y[k]-self.v0[k]*z[i,k] <= 1-x[i])
				Ext.addConstr(z[i,k] <= y[k])
				Ext.addConstr((self.v0[k]+self.v[i,k])*z[i,k] <= x[i])
		Ext.Params.timeLimit = tLimit
		Ext.update()
		Ext.optimize()

		wrtStr = str(Ext.objval)+"\t"+str(Ext.ObjBound)+"\t"+str(time.time()-t0)+"\t"+str(Ext.NodeCount)
		f = open("Results/"+str(self.name)+"_Ext.txt","a")
		f.write(wrtStr)
		f.close

	def BuildPrimal(self):
		self.PrimalMaster = Model('Master')
		for s in range(self.Nscen):
			self.theta[s] = self.PrimalMaster.addVar(lb=-GRB.INFINITY,name="theta"+str(s))
		for i in range(self.Nfv):
			self.x[i] = self.PrimalMaster.addVar(lb=0.0,ub=1.0,name="x"+str(i))
		self.PrimalMaster.addConstr(quicksum(self.x[i] for i in range(self.Nfv)) <= self.c*self.Nfv)
		self.PrimalMaster.setObjective(quicksum(self.pVec[s]*self.theta[s] for s in range(self.Nscen))+quicksum(self.cVec[j]*self.x[j] for j in range(self.Nfv)))
		self.PrimalMaster.setParam( 'OutputFlag', False )
		self.PrimalMaster.modelSense = GRB.MINIMIZE
		self.PrimalMaster.update()
		self.addBenders(init=True)

	def SolveBendersSub(self,scen_id,x_input):
		for i in self.Ck[scen_id]:
			self.scenrc[scen_id].getConstrByName("con:1-xi"+str(i)).RHS = 1-x_input[i]
			self.scenrc[scen_id].getConstrByName("con:xi"+str(i)).RHS = x_input[i]
		self.scenrc[scen_id].update()
		self.scenrc[scen_id].optimize()
		ObjValue = self.scenrc[scen_id].objval
		subg = {}
		for i in range(self.Nfv):
			if i in self.Ck[scen_id]:
				subg[i] = -self.scenrc[scen_id].getConstrByName("con:1-xi"+str(i)).pi+self.scenrc[scen_id].getConstrByName("con:xi"+str(i)).pi
			else:
				subg[i] = 0
		const = self.scenrc[scen_id].objval - sum(subg[i]*x_input[i] for i in range(self.Nfv))
		Benderscoef = []
		for i in range(self.Nfv):
			Benderscoef.append(subg[i])
		Benderscoef.append(const)
		Benderscoef = np.array([Benderscoef])
		self.BendersCuts[scen_id] = np.append(self.BendersCuts[scen_id], Benderscoef, axis=0)
		return ObjValue,const,subg

	def addGeneralSparse(self,s,x_input,I,shift=[0]):
		t0 = time.time()
		if len(shift) == 1:
			shift = [0]*(self.Nfv+1)

		I.sort()
		IsubList = []

		for card in range(min(len(I),int(self.c*self.Nfv))+1):
			IsubList = IsubList+list(itertools.combinations(I,card))

		SepProb = Model()
		mu = SepProb.addVars(I, lb=-GRB.INFINITY, name="mu")
		gamma = SepProb.addVar(lb=-GRB.INFINITY, name="gamma")
		SepProb.setObjective(quicksum(x_input[i]*mu[i] for i in I)+gamma)
		SepProb.modelSense = GRB.MAXIMIZE
		SepProb.setParam( 'OutputFlag', False )
		solns = {}
		Ibar = [i for i in range(self.Nfv) if i not in I]
		for i in Ibar:
			self.scensub[s].getVarByName("x"+str(i)).obj = -shift[i]
		for i in I:
			self.scensub[s].getVarByName("x"+str(i)).obj = 0
		for subI in IsubList:
			x_dict = {}
			for i in I:
				if i in subI:
					x_dict[i] = 1
				else:
					x_dict[i] = 0
			objx = 0.0
			if (tuple(I),subI) in self.ExactlyEvaluated[s].keys() and any(shift) == False:
				objx = self.ExactlyEvaluated[s][tuple(I),subI]
			else:
				rhsT0 = time.time()
				
				for i in I:
					if i in subI:
						self.scensub[s].getVarByName("x"+str(i)).lb = 1
					else:
						self.scensub[s].getVarByName("x"+str(i)).ub = 0
				self.scensub[s].update()
				self.scensub[s].optimize()
				objx = self.scensub[s].objval
				self.RhsT += time.time()-rhsT0
				
				if any(shift) == True and objx < sum(shift[i]*x_dict[i] for i in I)+shift[self.Nfv]:
					objx = sum(shift[i]*x_dict[i] for i in I)+shift[self.Nfv]
				
				if any(shift) == False:
					self.ExactlyEvaluated[s][tuple(I),subI] = objx
				for i in I:
					if i in subI:
						self.scensub[s].getVarByName("x"+str(i)).lb = 0
					else:
						self.scensub[s].getVarByName("x"+str(i)).ub = 1
				self.scensub[s].update()
			SepProb.addConstr(quicksum(x_dict[i]*mu[i] for i in I)+gamma<=objx)
		SepProb.update()
		SepProb.optimize()
		mu_soln = {}
		for i in I:
			mu_soln[i] = mu[i].x
		gamma_soln = gamma.x
		self.GeneralSparseT += time.time()-t0
		return SepProb.objval,mu_soln,gamma_soln,solns