from numpy import*
import matplotlib.pyplot as plt
import matplotlib

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)

plt.figure(1)
plt.loglog()
red    = [0.3,0.55]
dmv    = [-0.24,-0.33,-0.39,-0.48]
colors = ['blue','green','cyan','lime']
lbls   = ['R=70','R=55','R=70, R/3','R=55, R/3']
for z in range(len(red)):
	plt.figure(z)
	plt.loglog()
	for d in range(len(dmv)):
		data = loadtxt("M_dndlnM_ST_voids_Mpc-3_h0p6831_dmv"+str(dmv[d])+"_z"+str(red[z])+".txt").T
		plt.plot(data[0],data[1],ls='-',color=colors[d],label=r''+lbls[d])
	plt.xlim(1e2,1e16)
	plt.ylim(1e-25,5e6)
	plt.text(1e3,1e-15,r'z='+str(red[z]),fontsize=12)
	plt.xlabel(r'M [$M_\odot$]',fontsize=12)
	plt.ylabel(r'dn/dlnM [Mpc$^{-3}$]',fontsize=12)
	plt.legend(loc=0,prop={'size':12},frameon=False)


z03  = loadtxt("M_dndlnM_ST_voids_Mpc-3_h0p6831_dmv-0.24_z0.3.txt").T
z055 = loadtxt("M_dndlnM_ST_voids_Mpc-3_h0p6831_dmv-0.24_z0.55.txt").T
plt.figure(2)
plt.loglog()
plt.plot(z03[0],z03[1],'g--',label=r'z=0.3')
plt.plot(z055[0],z055[1],'b--',label=r'z=0.55')
plt.xlim(1e2,1e16)
plt.ylim(1e-25,5e6)
plt.xlabel(r'M [$M_\odot$]',fontsize=12)
plt.ylabel(r'dn/dlnM [Mpc$^{-3}$]',fontsize=12)
plt.legend(loc=0,prop={'size':12},frameon=False)
plt.show()  
