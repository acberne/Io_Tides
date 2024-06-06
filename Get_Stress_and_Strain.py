
import numpy as np


Love=np.loadtxt('output_Io.dat')

## Tidal Stresses on Enceladus
## Wahr et al (2009), Equations 12 - 19


def Get_Stress(thet,phi, M,h2=Love[4],l2=Love[6]):

	thet=thet*180/np.pi
	phi=phi*180/np.pi
	a_star = 421700e3 #Semi-major axis in m
	m_star = 1.8982e27 # Mass of Parent planet in kg
	G = 6.67e-11 #Gravitational constant
	Rs = 1821.6e3#Outer radius of the satellite in m
	kappa = 30e9 #Bulk modulus in Pa
	mu = 60e9 #Shear modulus in Pa
	lamb= kappa - 2*mu/3#lame parameter
	h_love= h2 # h Love number
	l_love= l2 # l Love number
	g = 1.796 # m/s^2
	n = 2*np.pi/(1.769*(24*3600)) ##Orbital angular velocity in s^-1
	eps = 0.0041 #Io oribital eccentricity

	t = (M*np.pi/180) #Time in radians

	alpha= (3*lamb+2*mu)/(lamb+2*mu)
	beta_1 = mu*(alpha*(h_love-3*l_love)+3*l_love)
	gamma_1= mu*(alpha*(h_love-3*l_love)-l_love)
	beta_2 = mu*(alpha*(h_love-3*l_love)-3*l_love)
	gamma_2= mu*(alpha*(h_love-3*l_love)+l_love)	

	Z = 3*G*m_star*Rs**2/(2*a_star**3)

	tau_thetthet = Z/(2*g*Rs) * (3*eps*(beta_1-gamma_1*np.cos(2*thet))*np.cos(t)*np.cos(2*phi)-eps*(beta_1+3*gamma_1*np.cos(2*thet))*np.cos(t)+4*eps*(beta_1-gamma_1*np.cos(2*thet))*np.sin(t)*np.sin(2*phi))

	tau_phiphi = Z/(2*g*Rs) * (3*eps*(beta_2-gamma_2*np.cos(2*thet))*np.cos(t)*np.cos(2*phi)-eps*(beta_2+3*gamma_2*np.cos(2*thet))*np.cos(t)+4*eps*(beta_2-gamma_2*np.cos(2*thet))*np.sin(t)*np.sin(2*phi))


	tau_phithet = 2*Z*l_love*mu/(g*Rs) *(4*eps*np.sin(t)*np.cos(thet)*np.cos(2*phi)-3*eps*np.cos(t)*np.cos(thet)*np.sin(2*phi)) 

	#Let's get strain tensor components as well.

	E = (9*kappa*mu)/(3*kappa+mu)
	nu = (3*kappa-2*mu)/(6*kappa+2*mu)
	eps_thetthet=(tau_thetthet-nu*tau_phithet)/E
	eps_phiphi=(-nu*tau_thetthet+tau_phithet)/E
	eps_thetphi=((1+nu)*tau_phithet)/E

	return tau_thetthet, tau_phiphi, tau_phithet, eps_thetthet, eps_phiphi, eps_thetphi



#print(Get_Stress(np.pi,np.pi/2,360)) #Arguments are Co-latitude in degrees, Longitude in degrees, Mean Anomaly in Degrees, h2, and l2, Output is Sigma_thetthet, Sigma_PhiPhi, Sigma_thetphi



