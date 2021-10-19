import master
import snipMaster
import capMaster
import llaMaster
import time

def SnipRoot(instanceName,Rbudget,snipNumber,K):
	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='cutpl')
	ss.AddDynamicSparse(K)
	ss.FreeMemory()

	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='cutpl')
	ss.GreedyViaSparse1(K)
	ss.FreeMemory()

	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='cutpl')
	ss.AddRandomSparse(K)
	ss.FreeMemory()

def CapRoot(instanceName,K,Nscen=250):
	ss = capMaster.CAPinst(instanceName,Nscen)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='level')
	ss.AddDynamicSparse(K)
	ss.FreeMemory()

	ss = capMaster.CAPinst(instanceName,Nscen)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='level')
	ss.GreedyViaSparse1(K)
	ss.FreeMemory()

	ss = capMaster.CAPinst(instanceName,Nscen)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='level')
	ss.AddRandomSparse(K)
	ss.FreeMemory()

def LlaRoot(Z,oneOverc,K,Cksize,kappa,n=500):
	ss = llaMaster.LLAinst(n,kappa,Z,oneOverc,Cksize)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='cutpl')
	ss.AddDynamicSparse(K)
	ss.FreeMemory()

	ss = llaMaster.LLAinst(n,kappa,Z,oneOverc,Cksize)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='cutpl')
	ss.GreedyViaSparse1(K)
	ss.FreeMemory()

	ss = llaMaster.LLAinst(n,kappa,Z,oneOverc,Cksize)
	ss.readData()
	ss.BuildPrimal()
	ss.addBenders(method='cutpl')
	ss.AddRandomSparse(K)
	ss.FreeMemory()

def SnipExtensive(instanceName,Rbudget,snipNumber):
	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.SolveExtensive()
	ss.FreeMemory()

def SnipBendersBnC(instanceName,Rbudget,snipNumber):
	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.BuildPrimal()
	tBenders = ss.addBenders(method='cutpl',store=True)
	ss.BendersBnC(methodName='Nonbasic',tLimit=60*60-tBenders)
	ss.FreeMemory()

	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.BuildPrimal()
	tBenders = ss.addBenders(method='cutpl',store=True)
	ss.BendersBnC(methodName='Basic',tLimit=60*60-tBenders)
	ss.FreeMemory()

def CapExtensive(instanceName,Nscen=250):
	ss = capMaster.CAPinst(instanceName,Nscen)
	ss.readData()
	ss.SolveExtensive()
	ss.FreeMemory()

def CapBendersBnC(instanceName,Nscen=250):
	ss = capMaster.CAPinst(instanceName,Nscen)
	ss.readData()
	ss.BuildPrimal()
	tBenders = ss.addBenders(method='level',store=True)
	ss.BendersBnC(methodName='Nonbasic',tLimit=60*60-tBenders)
	ss.FreeMemory()

def LlaExtensive(Z,oneOverc,Cksize,kappa,n=500):
	ss = llaMaster.LLAinst(n,kappa,Z,oneOverc,Cksize)
	ss.readData()
	ss.SolveExtensive()
	ss.FreeMemory()

def LlaBendersBnC(Z,oneOverc,Cksize,kappa,n=500):
	ss = llaMaster.LLAinst(n,kappa,Z,oneOverc,Cksize)
	ss.readData()
	ss.BuildPrimal()
	tBenders = ss.addBenders(method='cutpl',store=True)
	ss.BendersBnC(methodName='Nonbasic',tLimit=60*60-tBenders)
	ss.FreeMemory()

def SnipIsparseBnC(instanceName,Rbudget,snipNumber):
	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.BuildPrimal()
	t0 = time.time()
	ss.addBenders(method='cutpl')
	ss.AddDynamicSparse(4)
	tRoot = time.time()-t0
	ss.BendersBnC(methodName='IsparseNonbasic',tLimit=60*60-tRoot)
	ss.FreeMemory()

	ss = snipMaster.SNIPinst(instanceName,R=Rbudget,snipNo=snipNumber)
	ss.readData()
	ss.BuildPrimal()
	t0 = time.time()
	ss.addBenders(method='cutpl')
	ss.AddDynamicSparse(4)
	tRoot = time.time()-t0
	ss.BendersBnC(methodName='IsparseBasic',tLimit=60*60-tRoot)
	ss.FreeMemory()

def CapIsparseBnC(instanceName,Nscen=250):
	ss = capMaster.CAPinst(instanceName,Nscen)
	ss.readData()
	ss.BuildPrimal()
	t0 = time.time()
	ss.addBenders(method='level')
	ss.GreedyViaSparse1(7)
	tRoot = time.time()-t0
	ss.BendersBnC(methodName='IsparseNonbasic',tLimit=60*60-tRoot)
	ss.FreeMemory()

def LlaIsparseBnC(Z,oneOverc,Cksize,kappa,n=500):
	ss = llaMaster.LLAinst(n,kappa,Z,oneOverc,Cksize)
	ss.readData()
	ss.BuildPrimal()
	t0 = time.time()
	ss.addBenders(method='cutpl')
	ss.GreedyViaSparse1(7)
	tRoot = time.time()-t0
	ss.BendersBnC(methodName='IsparseNonbasic',tLimit=60*60-tRoot)
	ss.FreeMemory()

