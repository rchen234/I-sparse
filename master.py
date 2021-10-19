import numpy as np
import ast
from gurobipy import *
import time
import math
import itertools
import random

class StochIPinst:
	def __init__(self,instance):
		pass

	def addBenders(self,init=False,method='cutpl',tol=1e-4,timeLimit=60*60,store=False,ifBnC=False):
		if init == False:
			wrtStr = 'method='+str(method)+'\ttol='+str(tol)+'\ttimeLimit='+str(timeLimit)+'\n'
			if ifBnC == False:
				fileName = "Results/"+str(self.name)+"_Benders.txt"
			else:
				fileName = "Results/"+str(self.name)+"_BendersBnC.txt"
			f = open(fileName,"a")
			f.write(wrtStr)
			f.close

		t0 = time.time()
		TimeMasterLP = 0.0
		TimeCutLP = None
		TimeCutQP = None
		TimeBaseMIP = None
		TimeSub = None
		AvgBaseSize = None
		NscenSolved = None
		MaxTime = None
		MinTime = None
		MaxSize = None
		MinSize = None
		x_value = {}
		theta_value = {}
		if init == True:
			self.PrimalMaster.setObjective(0.0)
			self.PrimalMaster.update()
			self.PrimalMaster.optimize()
			for i in range(self.Nfv):	
				x_value[i] = self.x[i].x
			for s in range(self.Nscen):
				theta_value[s] = -float('inf')
			self.PrimalMaster.setObjective(quicksum(self.cVec[i]*self.x[i] for i in range(self.Nfv))+quicksum(self.pVec[s]*self.theta[s] for s in range(self.Nscen)))
			self.PrimalMaster.update()
		else:
			self.PrimalMaster.update()
			tStart = time.time()
			self.PrimalMaster.optimize()
			TimeMasterLP = TimeMasterLP+(time.time()-tStart)
			for i in range(self.Nfv):
				x_value[i] = self.x[i].x
			for s in range(self.Nscen):
				theta_value[s] = self.theta[s].x
		
		BendersUB = float('inf')
		LB = -float('inf')
		ContinueCondition = True
		iter = 0

		while ContinueCondition == True and time.time()-t0 < timeLimit:
			iter += 1
			NscenSolved = 0
			MaxTime = -float('inf')
			MinTime = float('inf')
			BendersCutsAdded = 0
			ContinueCondition = False
			CurObj = sum(self.cVec[i]*x_value[i] for i in range(self.Nfv))
			for s in range(self.Nscen):
				NscenSolved = NscenSolved+1
				tScen = time.time()
				ObjV,const,subg = self.SolveBendersSub(scen_id=s,x_input=x_value)
				if theta_value[s] < ObjV-tol*(abs(theta_value[s])+1) or init == True:
					self.thetaCutList[s].append(self.PrimalMaster.addConstr(self.theta[s] >= const+quicksum(subg[i]*self.x[i] for i in range(self.Nfv))))
					nzsubglist = [i for i in subg.keys() if subg[i]!=0]
					for i in nzsubglist:
						if i not in self.NonzeroInd[s]:
							self.NonzeroInd[s].append(i)
					abssub = [(abs(subg[i]),i) for i in subg.keys() if abs(subg[i])>=1e-6]
					abssub.sort(key = lambda x: x[0], reverse=True)
					indList = sorted([abssub[i][1] for i in range(min(6,len(abssub)))])
					if len(abssub) > 1 and indList not in self.NonzeroSparse[s]:
						self.NonzeroSparse[s].append(indList)

					self.Ncuts += 1
					BendersCutsAdded += 1
					cut = subg.copy()
					cut[self.Nfv] = const
					self.cutlist[s].append(cut)
					coef = subg.copy()
					coef[self.Nfv] = 1
					self.coeflist[s].append(coef)
					if method == 'cutpl':
						ContinueCondition = True
				CurObj = CurObj+self.pVec[s]*ObjV
				tScen = time.time()-tScen
				if tScen > MaxTime:
					MaxTime = tScen
				if tScen < MinTime:
					MinTime = tScen
			self.PrimalMaster.update()
			tStart = time.time()
			self.PrimalMaster.optimize()
			TimeMasterLP = TimeMasterLP+(time.time()-tStart)
			LB = self.PrimalMaster.objval
			for s in range(self.Nscen):
				theta_value[s] = self.theta[s].x
	
			if CurObj < BendersUB:
				BendersUB = CurObj

			if method == 'level' and BendersUB > LB+tol*(abs(LB)+1):
				if iter == 1 or sum(abs(x_value[i]-x_valueOld[i]) for i in range(self.Nfv)) > 1e-5:
					ContinueCondition = True

			if method == 'cutpl' and ContinueCondition == True:
				for j in range(self.Nfv):
					x_value[j] = self.x[j].x

			if method == 'level' and ContinueCondition == True:
				lt = LB+0.3*(BendersUB-LB)
				ltConstr = self.PrimalMaster.addConstr(quicksum(self.pVec[s]*self.theta[s] for s in range(self.Nscen))+quicksum(self.cVec[j]*self.x[j] for j in range(self.Nfv)) <= lt)
				self.PrimalMaster.setObjective(quicksum((self.x[j]-x_value[j])*(self.x[j]-x_value[j]) for j in range(self.Nfv)))
				self.PrimalMaster.update()
				self.PrimalMaster.optimize()
				x_valueOld = x_value.copy()
				for i in range(self.Nfv):
					x_value[i] = self.x[i].x
				self.PrimalMaster.remove(ltConstr)
				self.PrimalMaster.setObjective(quicksum(self.pVec[s]*self.theta[s] for s in range(self.Nscen))+quicksum(self.cVec[j]*self.x[j] for j in range(self.Nfv)))
				self.PrimalMaster.update()
				
			print('Benders Iter '+str(iter)+', PrimalMaster LB: '+str(LB)+', UB: '+str(BendersUB)+', Cuts Added: '+str(BendersCutsAdded))

			if init == True and iter >= 1:
				ContinueCondition = False

			if store == True:
				self.PrimalMaster.update()
				tStart = time.time()
				self.PrimalMaster.optimize()
				TimeMasterLP = TimeMasterLP+(time.time()-tStart)
				LB = self.PrimalMaster.objval
				wrtStr = str(iter)+'\t'+str(time.time()-t0)+'\t'+str(LB)+'\t'+str(TimeMasterLP)+'\t'+str(TimeCutLP)+'\t'+str(TimeCutQP)+'\t'+str(TimeBaseMIP)+\
					'\t'+str(TimeSub)+'\t'+str(AvgBaseSize)+'\t'+str(NscenSolved)+'\t'+str(MaxTime)+'\t'+str(MinTime)+'\t'+str(MaxSize)+'\t'+str(MinSize)+'\t'+\
					str(method)+'\t'+str(self.Ncuts)+'\t'+str(self.Nsubs)+'\n'
				f = open(fileName,"a")
				f.write(wrtStr)
				f.close
		print("LP time: "+str(time.time()-t0))
		return TimeMasterLP

	def IfromBenders(self,s,x_input,K,Exact=False,shift=[0]):
		I = []
		if any(shift) == False:
			shift = [0]*(self.Nfv+1)
		atleastK = False
		if Exact == False:
			theta_val = []
			theta_order = []
			for cut_ind in range(len(self.cutlist[s])):
				cut_order = 0
				cut = self.cutlist[s][cut_ind]
				theta_val.append(sum(cut[i]*x_input[i] for i in range(self.Nfv))+cut[self.Nfv])
				
			theta_order = list(np.argsort(theta_val))
			theta_order.reverse()
			find = 0
			while atleastK == False and find < len(self.cutlist[s]):
				ind = theta_order[find]
				find += 1
				nonzeros = [i for i in range(self.Nfv) if (self.cutlist[s][ind][i]-shift[i])*x_input[i]-min(self.cutlist[s][ind][i]-shift[i],0) > 1e-8]
				Ip = [i for i in range(self.Nfv) if i in I or i in nonzeros]
				if len(Ip) >= K:
					atleastK = True
					cut_max = self.cutlist[s][ind]
					inds = [i for i in nonzeros if i not in I]
				else:
					I = Ip.copy()
		else:
			atleastK = True
			ObjV,const,cut_max = self.SolveBendersSub(scen_id=s,x_input=x_input)
			inds = [i for i in range(self.Nfv) if cut_max[i] != 0]
		if atleastK == True:
			Delta = []
			for i in inds:
				Delta.append(-(cut_max[i]-shift[i])*x_input[i]+min(cut_max[i]-shift[i],0))
			Delta = np.array(Delta)
			sort_index = np.argsort(Delta)
			I = I+[inds[sort_index[k]] for k in range(K-len(I))]
			barI = [inds[sort_index[k]] for k in range(K-len(I),len(inds))]+[i for i in range(self.Nfv) if i not in inds and i not in I]
		else:
			barI = [i for i in range(self.Nfv) if i not in I]
		return sorted(I), barI
		

	def AddDynamicSparse(self,K,Exact=False,ifBendersApprox=0,ifshift=False,tLimit=30*60):
		t0 = time.time()		
		StopCondt = False
		IterNum = 0
		wrtStr = ""

		while StopCondt == False:
			StopCondt = True
			IterNum += 1
			x_soln = {}
			theta_soln = {}
			for i in range(self.Nfv):
				x_soln[i] = max(min(self.x[i].x,1),0)
			for s in range(self.Nscen):
				theta_soln[s] = self.theta[s].x
			for s in range(self.Nscen):
				if ifshift == True:
					bestcoef = None
					maxobj = -float("inf")
					for coef in self.BendersCuts[s]:
						if maxobj < sum(coef[i]*x_soln[i] for i in range(self.Nfv))+coef[self.Nfv]:
							maxobj = sum(coef[i]*x_soln[i] for i in range(self.Nfv))+coef[self.Nfv]
							bestcoef = coef
					
					for coef in self.TiltedIneqs[s]:
						if maxobj < sum(coef[i]*x_soln[i] for i in range(self.Nfv))+coef[self.Nfv]:
							maxobj = sum(coef[i]*x_soln[i] for i in range(self.Nfv))+coef[self.Nfv]
							bestcoef = coef
					
					shift = bestcoef
				else:
					shift = [0]*(self.Nfv+1)
				
				I, barI = self.IfromBenders(s,x_soln,K,Exact,shift)
				
				obj,mu_soln,gamma_soln,solns = self.addGeneralSparse(s,x_soln,I,shift=shift)

				newcoef = []
				for i in range(self.Nfv):
					if i in I:
						newcoef.append(mu_soln[i])
					else:
						newcoef.append(shift[i])
				newcoef.append(gamma_soln)
				Ibar = [i for i in range(self.Nfv) if i not in I]
				if obj+sum(shift[i]*x_soln[i] for i in Ibar) > self.theta[s].x+1e-6:
					self.TiltedIneqs[s].append(newcoef)
					if any(shift) == True:
						self.cutlist[s].append(newcoef)
					self.PrimalMaster.addConstr(self.theta[s] >= quicksum(newcoef[i]*self.x[i] for i in range(self.Nfv))+newcoef[self.Nfv])
					StopCondt = False
				if time.time()-t0 > tLimit:
					StopCondt = True
					break
			self.PrimalMaster.update()
			self.PrimalMaster.optimize()
			newStr = "ObjLB after iteration "+str(IterNum)+": "+str(self.PrimalMaster.objval)+", time: "+str(round(time.time()-t0))+", CGtime: "+str(round(self.GeneralSparseT))+", RHStime: "+str(round(self.RhsT))
			print(newStr)
			wrtStr = wrtStr+newStr+"\n"

		f = open("SparseResults/"+str(self.name)+"_"+str(K)+"_Dynamic.txt","a")
		f.write(wrtStr)
		f.close

	def GreedyViaSparse1(self,K,Exact=False,allcd=False,tLimit=30*60):
		t0 = time.time()
		shift = [0]*(self.Nfv+1)
		x_soln = {}
		theta_soln = {}
		for i in range(self.Nfv):
			x_soln[i] = max(min(self.x[i].x,1),0)
		Sparse1Cuts = {}
		LB = {}
		Flip = {}
		ArgSort = {}
		ScenCds = {}
		wrtStr = ""
		for s in range(self.Nscen):
			if allcd == False:
				ScenCds[s] = self.NonzeroInd[s]
			else:
				ScenCds[s] = range(self.Nfv)
			LB[s] = -float("inf")
			Sparse1Cuts[s] = {}
			for k in ScenCds[s]:
				I = [k]
				obj,mu_soln,gamma_soln,solns = self.addGeneralSparse(s,x_soln,I)
				if LB[s] < min(gamma_soln,mu_soln[k]+gamma_soln):
					LB[s] = min(gamma_soln,mu_soln[k]+gamma_soln)
				Sparse1Cuts[s][k] = (mu_soln[k],gamma_soln)
			for k in ScenCds[s]:
				obj0 = max(LB[s],Sparse1Cuts[s][k][1])
				obj1 = max(LB[s],Sparse1Cuts[s][k][0]+Sparse1Cuts[s][k][1])
				if obj0 <= obj1:
					Flip[s,k] = False
					minobj = obj0
					maxobj = obj1
				else:
					Flip[s,k] = True
					minobj = obj1
					maxobj = obj0
				Sparse1Cuts[s][k] = maxobj-minobj
			ArgSort[s] = sorted(Sparse1Cuts[s], key=Sparse1Cuts[s].get)
		# Sparsity 1 cuts not added yet
		print("Generate Sparsity 1 cuts: "+str(time.time()-t0))

		StopCondt = False
		IterNum = 0

		def maxViol(coef,x_soln,IndList):
			# Assume coef >= 0
			# Sort IndList by coef
			def takeCoef(elem):
				return coef[elem]
			IndList.sort(key=takeCoef)
			# Construct mapping: i -> index greater than i in IndList that has the largest x_soln[i] value
			NextLargest = {}
			NextLargest[IndList[-1]] = None
			ind = IndList[-1]
			maxVal = x_soln[ind]
			for k in range(2,len(IndList)+1):
				NextLargest[IndList[-k]] = ind
				if x_soln[IndList[-k]] > maxVal:
					ind = IndList[-k]
					maxVal = x_soln[IndList[-k]]
			I = []
			I.append(ind)
			while NextLargest[ind] != None:
				ind = NextLargest[ind]
				I.append(ind)
			vio = coef[I[0]]*x_soln[I[0]]
			Iind = 1
			while Iind < len(I):
				vio += (coef[I[Iind]]-coef[I[Iind-1]])*x_soln[I[Iind]]
				Iind += 1
			return vio

		while StopCondt == False:
			StopCondt = True
			IterNum += 1
			for i in range(self.Nfv):
				x_soln[i] = max(min(self.x[i].x,1),0)
			for s in range(self.Nscen):
				theta_soln[s] = self.theta[s].x
			for s in range(self.Nscen):
				I = []
				x_flipped = {k:x_soln[k]+Flip[s,k]*(1-2*x_soln[k]) for k in ScenCds[s]}
				while len(I) < min(K,len(ScenCds[s])):
					Ibar = [i for i in ScenCds[s] if i not in I]
					iBest = None
					vio_max = -float("inf")
					for i in Ibar:
						Icand = I.copy()
						Icand.append(i)
						IcandViol = maxViol(Sparse1Cuts[s],x_flipped,Icand)
						if IcandViol > vio_max:
							vio_max = IcandViol
							iBest = i
					I.append(iBest)

				obj,mu_soln,gamma_soln,solns = self.addGeneralSparse(s,x_soln,I)

				newcoef = []
				for i in range(self.Nfv):
					if i in I:
						newcoef.append(mu_soln[i])
					else:
						newcoef.append(shift[i])
				newcoef.append(gamma_soln)
				Ibar = [i for i in range(self.Nfv) if i not in I]
				if obj+sum(shift[i]*x_soln[i] for i in Ibar) > self.theta[s].x+1e-6:
					self.TiltedIneqs[s].append(newcoef)
					if any(shift) == True:
						self.cutlist[s].append(newcoef)
					self.PrimalMaster.addConstr(self.theta[s] >= quicksum(newcoef[i]*self.x[i] for i in range(self.Nfv))+newcoef[self.Nfv])
					StopCondt = False
				if time.time()-t0 > tLimit:
					StopCondt = True
					break
			self.PrimalMaster.update()
			self.PrimalMaster.optimize()
			newStr="ObjLB after iteration "+str(IterNum)+": "+str(self.PrimalMaster.objval)+", time: "+str(round(time.time()-t0))+", CGtime: "+str(round(self.GeneralSparseT))+", RHStime: "+str(round(self.RhsT))
			print(newStr)
			wrtStr = wrtStr+newStr+"\n"

		f = open("SparseResults/"+str(self.name)+"_"+str(K)+"_Greedy.txt","a")
		f.write(wrtStr)
		f.close

	def AddRandomSparse(self,K,tLimit = 30*60):
		t0 = time.time()
		StopCondt = False
		IterNum = 0
		random.seed(0)
		wrtStr = ""
		while StopCondt == False:
			IterNum += 1
			x_soln = {}
			theta_soln = {}
			for i in range(self.Nfv):
				x_soln[i] = max(min(self.x[i].x,1),0)
			for s in range(self.Nscen):
				theta_soln[s] = self.theta[s].x
			for s in range(self.Nscen):
				I = random.sample(self.NonzeroInd[s], min(K,len(self.NonzeroInd[s])))
				obj,mu_soln,gamma_soln,solns = self.addGeneralSparse(s,x_soln,I)
				if obj > self.theta[s].x+1e-6:
					self.PrimalMaster.addConstr(self.theta[s] >= quicksum(mu_soln[i]*self.x[i] for i in I)+gamma_soln)
					StopCondt = False
				if time.time()-t0 > tLimit:
					StopCondt = True
					break
			self.PrimalMaster.update()
			self.PrimalMaster.optimize()
			newStr = "ObjLB after iteration "+str(IterNum)+": "+str(self.PrimalMaster.objval)+", time: "+str(round(time.time()-t0))+", CGtime: "+str(round(self.GeneralSparseT))+", RHStime: "+str(round(self.RhsT))
			print(newStr)
			wrtStr = wrtStr+newStr+"\n"

		f = open("SparseResults/"+str(self.name)+"_"+str(K)+"_Random.txt","a")
		f.write(wrtStr)
		f.close

	def BendersBnC(self,methodName,tol=1e-6,tLimit=60*60):
		t0 = time.time()
		global tBenders
		global MIPSOLcount
		tBenders = 0.0
		MIPSOLcount = 0
		for i in range(self.Nfv):
			self.x[i].vtype = GRB.BINARY
		self.PrimalMaster.setParam( 'OutputFlag', True )
		self.PrimalMaster.Params.LogFile="Results/"+str(self.name)+"_"+methodName+"_BendersBnC.log"
		self.PrimalMaster.update()

		if methodName == "Basic" or "IsparseBasic":
			self.PrimalMaster.Params.Presolve = 0
			self.PrimalMaster.Params.Cuts = 0

		def cb(model, where):
			if where == GRB.Callback.MIPSOL:
				global tBenders
				global MIPSOLcount
				MIPSOLcount = MIPSOLcount+1
				x_value = {}
				theta_value = {}
				for i in range(self.Nfv):
					x_value[i] = max(min(model.cbGetSolution(model.getVarByName("x"+str(i))),1),0)
				for s in range(self.Nscen):
					theta_value[s] = model.cbGetSolution(model.getVarByName("theta"+str(s)))
				for s in range(self.Nscen):
					tBendersStart = time.time()
					ObjV,const,subg = self.SolveBendersSub(scen_id=s,x_input=x_value)
					tBenders = tBenders+(time.time()-tBendersStart)
					if theta_value[s] < ObjV-tol*(abs(theta_value[s])+1):
						model.cbLazy(model.getVarByName("theta"+str(s)) >= const+quicksum(subg[i]*model.getVarByName("x"+str(i)) for i in range(self.Nfv)))
		self.PrimalMaster.Params.lazyConstraints = 1
		self.PrimalMaster.Params.TimeLimit = tLimit

		self.PrimalMaster.update()

		self.PrimalMaster.optimize(cb)
		print("Optimal Value: "+str(self.PrimalMaster.objval)+", B&C time: "+str(time.time()-t0))
		f = open("Results/"+str(self.name)+"_"+methodName+"_BendersBnC.txt","w")
		wrtStr = str(self.PrimalMaster.objval)+"\t"+str(self.PrimalMaster.ObjBound)+"\t"+str(time.time()-t0)+"\t"+str(self.PrimalMaster.NodeCount)\
			+"\t"+str(60*60-tLimit)+"\t"+str(self.GeneralSparseT)+"\t"+str(self.RhsT)
		f.write(wrtStr)
		f.close()