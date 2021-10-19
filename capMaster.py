import numpy as np
import ast
from gurobipy import *
import master
import time
import itertools
#from sympy.combinatorics.graycode import GrayCode

class CAPinst(master.StochIPinst):
	def __init__(self,instance,Nscen):
		self.Nscen = Nscen
		self.Nfv = None
		self.instance = instance
		self.Jcard = None
		self.Icard = None
		self.name = instance
		self.cVec = None # f vector
		self.pVec = None
		self.dVec = None
		self.TMtr = None
		self.WMtr = None
		self.hVec = None
		self.sVec = None
		self.b = None
		self.BendersCuts = {}
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
		self.scenrc = {}
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
		# Graph edges
		self.Edges = {}

	def _readVec(self, vec, nType=np.float32):
		return np.asarray(ast.literal_eval(vec[:-1]), dtype = nType)

	def _readMtr(self, mtr):
		mtr_list = []
		for line in mtr:
			mtr_list.append([float(y) for y in line.strip('\t\n[],').split(',')])
		return np.array(mtr_list)


	def readData(self):
		f = open('Instances/CAP/'+self.instance+'.dat','r')
		x = f.readlines()
		f.close()
		self.Jcard = 50
		if len(x) > 5200:
			self.Icard = 50
		else:
			self.Icard = 25
		self.Nfv = self.Icard
		self.pVec = [1/self.Nscen]*self.Nscen
		
		self.cVec = self._readVec(x[0], np.float32)
		self.dVec = self._readVec(x[1], np.float32)
		self.TMtr = self._readMtr(x[2:2+self.Icard+self.Jcard])
		self.WMtr = self._readMtr(x[2+self.Icard+self.Jcard:2+2*self.Icard+2*self.Jcard])
		self.hVec = self._readMtr(x[2+2*self.Icard+2*self.Jcard:self.Nscen+2+2*self.Icard+2*self.Jcard])
		self.sVec = self._readMtr([x[-1]])[0]
		sumlam = []
		for s in range(self.Nscen):
			sumlam.append(sum(self.hVec[s]))
		self.b = np.ceil(max(sumlam)/self.sVec[0])
		for s in range(self.Nscen):
			self.cutlist[s] = []
			self.coeflist[s] = []
			self.thetaCutList[s] = []
			self.BendersCuts[s] = np.empty((0,self.Nfv+1), float)
			self.NonzeroInd[s] = []
			self.NonzeroSparse[s] = []
			self.Solns[s] = []
			self.ExactlyEvaluated[s] = {}
			self.Edges[s] = {}
			self.IEvaluated[s] = []
			self.Iineqs[s] = []
			self.TiltedIneqs[s] = []
		self.Initialize()

	def Initialize(self):
		t0 = time.time()
		self.scenrc = Model()
		self.scenrc.setParam('OutputFlag', False)
		self.scenrc.modelSense = GRB.MINIMIZE
		y = self.scenrc.addVars(range(self.Icard*self.Jcard),name="y")
		scencon = self.scenrc.addConstrs((quicksum(self.WMtr[k][i]*y[i] for i in range(self.Icard*self.Jcard)) >= 0.0 for k in range(self.Jcard+self.Icard)),name="scencon")
		self.scenrc.setObjective(quicksum(self.dVec[i]*y[i] for i in range(self.Icard*self.Jcard)))
		self.scenrc.update()
		print("Build recourse subproblems: "+str(time.time()-t0))
		self.scensub = Model()
		self.scensub.setParam('OutputFlag', False)
		self.scensub.modelSense = GRB.MINIMIZE
		x = self.scensub.addVars(range(self.Nfv),lb=0.0,ub=1.0,name="x")
		y = self.scensub.addVars(range(self.Icard*self.Jcard),obj=self.dVec,name="y")
		self.scensubCon = self.scensub.addConstrs((quicksum(self.WMtr[k][i]*y[i] for i in range(self.Icard*self.Jcard))+quicksum(self.TMtr[k][i]*x[i] for i in range(self.Nfv)) >= 0.0 for k in range(self.Jcard+self.Icard)))
		#self.scensub.setObjective(quicksum(self.dVec[i]*y[i] for i in range(self.Icard*self.Jcard)))
		self.scensub.update()
		print("Build scenario subproblems: "+str(time.time()-t0))


	def BuildPrimal(self):
		for s in range(self.Nscen):
			self.theta[s] = self.PrimalMaster.addVar(lb=-GRB.INFINITY,name="theta"+str(s))
		for i in range(self.Nfv):
			self.x[i] = self.PrimalMaster.addVar(lb=0.0,ub=1.0,name="x"+str(i))
		self.PrimalMaster.addConstr(quicksum(self.x[i] for i in range(self.Nfv)) >= self.b)
		self.PrimalMaster.setObjective(quicksum(self.pVec[s]*self.theta[s] for s in range(self.Nscen))+quicksum(self.cVec[j]*self.x[j] for j in range(self.Nfv)))
		self.PrimalMaster.setParam( 'OutputFlag', False )
		self.PrimalMaster.modelSense = GRB.MINIMIZE
		self.PrimalMaster.update()
		self.addBenders(init=True)

	def FreeMemory(self):
		del self.scensub
		del self.scenrc
		del self.PrimalMaster

	def SolveBendersSub(self,scen_id,x_input,store=False):
		x_input = np.array([x_input[i] for i in range(self.Nfv)])
		RHS = np.dot(self.TMtr,x_input)
		for k in range(self.Jcard+self.Icard):
			self.scenrc.getConstrByName('scencon['+str(k)+']').rhs = -RHS[k]+self.hVec[scen_id][k]
		self.scenrc.update()
		self.scenrc.optimize()
		ObjValue = self.scenrc.objval
		dual_soln = []
		for k in range(self.Jcard+self.Icard):
			dual_soln.append(self.scenrc.getConstrByName('scencon['+str(k)+']').pi)
		dual_soln = np.array(dual_soln)
		prod = np.dot(dual_soln,self.TMtr)
		subg = {}
		Benderscoef = []
		for i in range(self.Nfv):
			subg[i] = -prod[i]
			Benderscoef.append(-prod[i])
		const = self.scenrc.objval - sum(subg[i]*x_input[i] for i in range(self.Nfv))
		Benderscoef.append(const)
		if store == True:
			soln = list(x_input.copy())
			soln.append(ObjValue)
			self.Solns[scen_id].append(soln)
		Benderscoef = np.array([Benderscoef])
		self.BendersCuts[scen_id] = np.append(self.BendersCuts[scen_id], Benderscoef, axis=0)
		return ObjValue,const,subg

	def addGeneralSparse(self,s,x_input,I,shift=[0]):
		t0 = time.time()
		if len(shift) == 1:
			shift = [0]*(self.Nfv+1)
		I.sort()
		IsubList = []
		for card in range(max(0,len(I)+self.b-self.Nfv),len(I)+1):
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
			self.scensub.getVarByName("x["+str(i)+"]").obj = -shift[i]
		for i in I:
			self.scensub.getVarByName("x["+str(i)+"]").obj = 0
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
				
				for k in range(self.Jcard+self.Icard):
					self.scensubCon[k].RHS = self.hVec[s][k]
				for i in I:
					if i in subI:
						self.scensub.getVarByName("x["+str(i)+"]").lb = 1
					else:
						self.scensub.getVarByName("x["+str(i)+"]").ub = 0
				self.scensub.update()
				self.scensub.optimize()
				objx = self.scensub.objval
				self.RhsT += time.time()-rhsT0
				
				if any(shift) == True and objx < sum(shift[i]*x_dict[i] for i in I)+shift[self.Nfv]:
					objx = sum(shift[i]*x_dict[i] for i in I)+shift[self.Nfv]
				
				if any(shift) == False:
					self.ExactlyEvaluated[s][tuple(I),subI] = objx
				for i in I:
					if i in subI:
						self.scensub.getVarByName("x["+str(i)+"]").lb = 0
					else:
						self.scensub.getVarByName("x["+str(i)+"]").ub = 1
				self.scensub.update()
			SepProb.addConstr(quicksum(x_dict[i]*mu[i] for i in I)+gamma<=objx)
		SepProb.update()
		SepProb.optimize()
		mu_soln = {}
		for i in I:
			mu_soln[i] = mu[i].x
		gamma_soln = gamma.x
		self.GeneralSparseT += time.time()-t0
		return SepProb.objval,mu_soln,gamma_soln,solns

	def SolveExtensive(self,tLimit=60*60):
		t0 = time.time()
		Ext = Model()
		Ext.Params.LogFile="Results/"+str(self.name)+"_Ext.log"
		x = Ext.addVars(self.Nfv, vtype=GRB.BINARY)
		y = Ext.addVars(self.Nscen,self.Icard,self.Jcard)
		Ext.addConstr(quicksum(x[i] for i in range(self.Nfv)) >= self.b)
		Ext.addConstrs(quicksum(y[s,i,j] for i in range(self.Icard)) >= self.hVec[s][j] for s in range(self.Nscen) for j in range(self.Jcard))
		Ext.addConstrs(-quicksum(y[s,i,j] for j in range(self.Jcard))+self.TMtr[self.Jcard+i,i]*x[i] >= 0 for s in range(self.Nscen) for i in range(self.Icard))
		Ext.setObjective(quicksum(self.cVec[i]*x[i] for i in range(self.Nfv))+quicksum(self.pVec[s]*self.dVec[i*self.Jcard+j]*y[s,i,j] for s in range(self.Nscen) for i in range(self.Icard) for j in range(self.Jcard)))
		Ext.Params.timeLimit = tLimit
		Ext.update()
		Ext.optimize()

		wrtStr = str(Ext.objval)+"\t"+str(Ext.ObjBound)+"\t"+str(time.time()-t0)+"\t"+str(Ext.NodeCount)
		f = open("Results/"+str(self.name)+"_Ext.txt","a")
		f.write(wrtStr)
		f.close