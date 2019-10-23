import numpy as np

def event(r,s,e,i):
	#expose
	if r == 0:
		s = s - 1.0
		e = e + 1.0
	#expose to infect
	elif r == 1:
		e = e - 1.0
		i = i + 1.0
	#recovery
	else:
		i = i - 1.0
	return s, e, i

def gillespie(R0,gamma,sigma,N,n_sim,t_end = 12*24.0):

	# use R0 to set infection scale, no births or immigration
	beta = R0*gamma
	t_lin = np.linspace (0,t_end,2000)

	m_interp = np.zeros((n_sim,t_lin.shape[0]))
	m_new_interp = np.zeros((n_sim,t_lin.shape[0]))
	n_interp = np.zeros((n_sim,t_lin.shape[0]))
	bins = np.linspace(0,t_end,t_end*gamma)
	m_new_bins = np.zeros((n_sim,bins.shape[0]))


	for ii in np.arange(n_sim):
		
		#initial number of infected
		e = 1
		m = 0
		#initial number of susceptibles
		n = N-m

		m_ts = []
		m_new_ts = []
		e_ts = []
		n_ts = []
		t_ts = []
		t_date = []

		t = 0.0

		while t < t_end:

		#update rates
			r_i = beta*n/N*m
			r_e = sigma*e
			r_r = gamma*m
			r = r_i + r_r + r_e

			if r < 1e-9:
				break
			#draw random numbers
			dt = -1.0/r*np.log(np.random.uniform(0,1))
			P = r*np.random.uniform(0,1)

			r_array = np.array([[r_i,0],[r_e,1],[r_r,2]])
			r_array = r_array[np.argsort(r_array[:, 0])]

			for jj in np.arange(r_array.shape[0]):
				if P < np.sum(r_array[0:jj+1,0]):
					n, e, m = event(r_array[jj,1],n,e,m)
					break
				
			if t > 0:
				if m-m_ts[-1]>0:
					m_new_ts.append(m-m_ts[-1])
				else:
					m_new_ts.append(0)
			else:
				m_new_ts.append(0)		
			
			m_ts.append(m)
			e_ts.append(e)
			n_ts.append(n)
			t_ts.append(t)

			#update times
			t = t + dt

			if m < 0.0 or n < 0.0:
				break

			for bb in np.arange(bins.shape[0]-1):
				if t < bins[bb+1] and t >  bins[bb]:
					m_new_bins[ii,bb] += m_new_ts[-1]

			if t > bins[-1]:
				m_new_bins[ii,-1] += m_new_ts[-1]

		m_interp[ii,:] = np.interp(t_lin,t_ts,m_ts)
		m_new_interp[ii,:] = np.interp(t_lin,t_ts,m_new_ts)
		n_interp[ii,:] = np.interp(t_lin,t_ts,n_ts)

	return m_new_bins