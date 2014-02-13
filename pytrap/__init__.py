#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division

import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import math
from numpy import log, exp, pi, sqrt, cosh, vectorize
import numpy as np
from scipy.integrate import odeint
from multiprocessing import Pool
from pylab import figure, imshow, xlabel, ylabel, plt, show, contour, meshgrid, cm, savefig, plot, xlim, ylim, scatter
import time
import os


# Global trap parameters
Rz=0.013 #radius of the einzel electrode
R=0.008 # radii of the other electodes
L = 0.22 # half the length of trap for normalization
omega=1.3152/pi #constant (see Bertram'method)
eps=0.0001 #for the differentiations
#positions of the electrodes in the W-space (calculated in mathematica see article)
xiw = np.array([10.9948, 14.1364, 15.9821, 17.5529, 18.5739, 21.2836, 22.4618, 24.6229, 25.8437, 28.4623, 30,
				31.9182, 34.4907, 36.4722])

# Shortcuts
xiw_left = xiw[:-1]
xiw_right = xiw[1:]
xiw_deltas = xiw[1:]-xiw[:-1]
xiw_sums = xiw[1:]+xiw[:-1]
bb = (Rz/R)**2

def ff(zp):
	"""
	function that goes from the Z-plane to the W-plane. (from mathematica)
	"""
	z = -abs(zp*L)+0.1525
	if z>= 0.007:
		return log((bb-1)/4*(exp(pi*z/Rz-1/sqrt(bb)*log((sqrt(bb)-1)/(sqrt(bb)+1)))+(1+bb)/(1-bb))) + 30
	elif -0.007<z<0.007:
		return -0.655414 + 321.236*z - 4238.41*z**2 + 512787.0*z**3 + (4.11932*10**6)*z**4 - (4.34084*10**9)*z**5 + 30
	else:
		return -log((bb-1)/(4*bb)*exp(-sqrt(bb)*(pi*z/Rz-log((sqrt(bb)+1)/(sqrt(bb)-1))))-(1+bb)/(4*bb)) + 30

class Setup(object):
		"""
		Defines the setup of the experiment.
		
		V -- list of electrode potentials. eg: [6500, 5850, 4150, 1650, 4950]
		ma, ne, ener -- atomic mass (ua), charge number and energy (eV) of the particles.
		alpha -- space charge coeficient
		b -- friction coefficient
		
						eg:     m = 16, q = 4,  ener = 5200.
		"""
		def __init__(self, V = [8000, 5850, 4150, 1650, 3300], ma = 16, ne = 4, ener = 5200, alpha=0, b = 0):
			self.V = V
			self.m = ma*1.672623e-27
			self.q = ne*1.602e-19
			self.ener = ener
			self.vz = sqrt(2*self.q*ener/self.m)
			self.T = L/self.vz #characteristic timescale used for adimentionning
			self.alpha = alpha
			self.b = b

class Trap(object):

	def __init__(self, setup):
		"""
		Creates a trap with a given setup.
		"""
		V = setup.V
		Vi = np.array([0, V[0], V[0], V[1], V[1], V[2], V[2], V[3], V[3], 0, 0, V[4], V[4], 0])
		
		self.V = setup.V
		self.ener = setup.ener
		# Optimizations:
		self.Vi_deltas = Vi[1:]-Vi[:-1]
		self.Vi_left = Vi[:-1]
		self.deltas_ratio = self.Vi_deltas/xiw_deltas

		# Description of the particles:
		self.q_over_m = setup.q/setup.m
		self.vz = setup.vz
		
		self.T = setup.T
		self.mass = setup.m
		self.eta = 1/(2*setup.ener)
		self.alpha = setup.alpha
		self.b = setup.b
		self.ener = setup.ener

	def Vw(self, w):
		"""
		w -- position on the axis in the W-plane
		returns the potential on the axis
		"""
		potentials = ((self.deltas_ratio*xiw_left-self.Vi_left)*np.cosh(omega*xiw_deltas)/(np.cosh(omega*(2*w+xiw_sums))+np.cosh(omega*xiw_sums))
			      -1/(2*omega)*self.deltas_ratio*np.log(np.cosh(omega*(w-xiw_right))/np.cosh(omega*(w-xiw_left))))
		return potentials.sum()

	def d2V(self, z): 
		"""
		calculates the second derivative of the potential along axis
		z has no dimension
		"""
		a = self.Vw(ff(z+2*eps))
		b = self.Vw(ff(z-2*eps))
		c = self.Vw(ff(z))
		return (a+b-2*c)/(4*eps**2)

	def potential(self, z, r):
		"""
		z, r -- position in trap (in meters)
		returns the electrostatic potential (order 2 of the development)
		"""
		return self.Vw(ff(z/L)) - r**2 * self.d2V(z)/4 

	def Ez0(self, z):
		"""
		calculates the electric field ON THE AXIS along z
		z has no dimension
		"""
		a=self.Vw(ff(z+eps))
		b=self.Vw(ff(z-eps))
		return -(a-b)/(2*eps)

	def fz(self):
		"""
		calcultates the frequency od the longitudinal oscillation movement
		"""
		time = np.linspace(0, 5, 1000)
		def sm1(x, t):
			if x[0] < 0:
				return np.array([x[1], 0]) 
			else:
				return np.array([x[1], self.eta*self.Ez0(x[0])])
		traj = odeint(sm1,[0,1.],time)
		if traj[-1,0] < 0:
			ind = np.where(traj[:,0] > 0)[0][-1]
			t = [time[ind+1],time[ind]]
			z = [traj[ind+1,0],traj[ind,0]] #must be increasing for interp
			return 1/(np.interp(0,z,t)*self.T)
		else:
			return -1

	def delta_r(self):
		"""
		returns the coeficient delta indicating the RADIAL stability of the movement
		"""
		time = np.linspace(0, 5, 1000)
		def sm1(x, t):
			if x[0] < 0:
				return np.array([x[1], 0]) 
			else:
				return np.array([x[1], self.eta*self.Ez0(x[0])])
		traj = odeint(sm1,[0,1.],time)
		if traj[-1,0] < 0:
			ind = np.where(traj[:,0] > 0)[0][-1]
			tm = [time[ind+1],time[ind]]
			zm = [traj[ind+1,0],traj[ind,0]]
			period = np.interp(0,zm,tm)
			def sm2(x, t):
				t = t - np.floor(t/period)*period
				pos = np.interp(t, time, traj[:,0])
				psi = self.alpha + self.b**2/4. + (self.eta/2.)*self.d2V(pos)
				return np.array([x[1], psi*x[0], x[3], psi*x[2]])
			x0 = [1, 0, 0, 1]
			sol = odeint(sm2, x0, np.linspace(0, period, 1000))
			delta = np.abs((sol[-1,0] + sol[-1,3])/2)
			beta  = np.arccos((sol[-1,0] + sol[-1,3])/2)/np.pi
			if self.b == 0.:
				return delta
			else:
				out = np.log(delta + np.sqrt(np.abs(delta**2-1)))/period - self.b
				return out
		else:
			return 1.0
			
	def delta_s(self):
		"""
		returns the coeficient delta indicating the LONGITUDINAL stability of the movement 
		"""
		time = np.linspace(0, 2.5, 1000)
		def sm1(x, t):
			return np.array([x[1], self.eta*self.Ez0(x[0])])
		traj = odeint(sm1,[0,1.],time)
		if traj[-1,0] < 0:
			ind = np.where(traj[:,0] > 0)[0][-1]
			tm = [time[ind+1],time[ind]]
			zm = [traj[ind+1,0],traj[ind,0]]
			period = np.interp(0,zm,tm)
			def sm2(x, t):
				t = t - np.floor(t/period)*period
				pos = np.interp(t, time, traj[:,0])
				psi = self.alpha - self.eta*self.d2V(pos)
				return np.array([x[1], psi*x[0], x[3], psi*x[2]])
			x0 = [1, 0, 0, 1]
			sol = odeint(sm2, x0, np.linspace(0, period, 1000))
			delta = np.abs((sol[-1,0] + sol[-1,3])/2)
			return delta
		else:
			return 1.0
			
	def poincare(self, num_turn, rmin, rmax, n_rad):
		"""
		Plot the Poincaré section of the given setup.
			num_turns: number of crossings of the Poincaré map.
			rmin, rmax: initial radius of the particles will be linearly distributed between
			n_rad: number of particles 
		"""
		potentials = 0.001*np.array(self.V)
		str_pot = ""
		str_ener = str(self.ener/1000.)
		for potential in potentials:
			str_pot += str(potential)+" "
		radii = np.linspace(rmin, rmax, n_rad)
		for radius in radii:
			cmd = "cytrap/solveode -n%s %s %s 0 %s >> out.dat" % (num_turn, str_pot, str_ener, radius)
			os.popen(cmd)
		data = np.loadtxt('out.dat') 
		t0, x0, y0, z0, px0, py0, pz0 = data.T
		os.popen('rm out.dat')
		r0 = np.sqrt(x0**2+y0**2)
		pr0 = (x0*px0+y0*py0)*r0**(-1)
		figure()
		xlim([-0.005, 0.005])
		ylim([-0.065,0.065])
		ylabel("$p_r$",fontsize="x-large")
		xlabel("$r$",fontsize="x-large")
		scatter(x0, pr0, s=0.3, color='k')
		show()
		
	def trajectory(self, tend, rmin, rmax, n_rad):
		"""
		Compute the trajectory of an ion. 
			tend: time of the end of the trajectory
			rmin, rmax: initial radius of the particles will be linearly distributed between
			n_rad: number of particles
		"""
		potentials = 0.001*np.array(self.V)
		str_pot = ""
		str_ener = str(self.ener/1000.)
		for potential in potentials:
			str_pot += str(potential)+" "
		radii = np.linspace(rmin, rmax, n_rad)
		figure()
		xlim([-0.22,0.22])
		ylim([-0.008,0.008])
		xlabel("z (m)")
		ylabel("x (m)")
		for radius in radii:
			cmd = "cytrap/solveode -t%s %s %s 0 %s >> out.dat" % (tend, str_pot, str_ener, radius)
			os.popen(cmd)
			data = np.loadtxt('out.dat') 
			t0, x0, y0, z0, px0, py0, pz0 = data.T
			os.popen('rm out.dat')
			plot(z0, x0)
		show()
			
class Zgraph(object):

	def __init__(self, setup, (V1s, V1e), (Vzs, Vze), nb = 8, mode = 'stab'):
		"""
		Creates the stability map for a given Setup.
		   __init__(self, setup, (V1s, V1e), (Vzs, Vze), nb = 8, mode = 'stab')

		setup: the setup
		(V1s, V1e): start and end points of the map along V1
		(Vzs, Vze): start and end points of the map along Vz
		nb: number of points along each side of the map
		mode: 'stab' gives the radial stability, 'sync' gives the synchronization stability
		"""
		self.setup = setup
		(self.V1s, self.V1e) = (V1s, V1e)
		(self.Vzs, self.Vze) = (Vzs, Vze)
		self.nb = nb
		self.VZ = np.linspace(Vzs, Vze, nb)
		self.V1 = np.linspace(V1s,V1e, nb)
		self.X, self.Y = meshgrid(self.VZ, self.V1)
		self.mode = mode
	
	def calc(self):
		p = Pool()
		TASKS = [(worker, (i, self.setup, self.X, self.Y, self.mode)) for i in range(self.nb)]
		self.Z = p.map(calculatestar, TASKS)
	
	def plot(self):
		Zdelta = []
		for aa in self.Z:
			Zdelta += [aa-1]
		figure()
		imshow(Zdelta, interpolation='bilinear', origin='lower',
	                	cmap=cm.bone, extent=(self.Vzs, self.Vze, self.V1s, self.V1e))
		CS = contour(self.X, self.Y, Zdelta, [0], linewidths=4, colors='white')
		CS2 = contour(self.X, self.Y, Zdelta, 16, linewidths=1, colors='k')
		plt.clabel(CS2, fontsize=6, inline=1)
		plt.clabel(CS, fontsize=9, inline=1)
		if self.mode == 'stab':
			plt.title('Radial stability map %s V\n %s\n alpha=%s\n b=%s' % (str(self.setup.ener), str(self.setup.V), str(self.setup.alpha), str(self.setup.b)))
		else:
			plt.title('Synchronization stability map %s V\n %s\n alpha=%s\n b=%s' % (str(self.setup.ener), str(self.setup.V), str(self.setup.alpha), str(self.setup.b)))
		xlabel('Vz (V)')
		ylabel('V1 (V)')
		show()

# the folowing functions are separted from the classes because they are called by multiprocessing
def worker(i, setup, X, Y, mode):
	line = stability(setup, Y[i], X[i], mode)
	return line
def calculate(func, args):
	return func(*args)	
def calculatestar(args):
	return calculate(*args)
@vectorize
def stability(setup, x, y, mode):
	setup.V[0] = x
	setup.V[4] = y
	trap_conf = Trap(setup)
	if mode == 'stab':
		d = trap_conf.delta_r()
	else:
		d = trap_conf.delta_s()
	del(trap_conf)
	return d
	
class movie(object):
	def __init__(self, setup, (V1s, V1e), (Vzs, Vze), nb, (astart, aend), (bstart, bend), nbf, path=""):
		"""
		Creates a animation of stability map.

		setup: the setup
		(V1s, V1e): start and end points of the map along V1
		(Vzs, Vze): start and end points of the map along Vz
		nb: number of points along each side of the map
		(astart, aend): starting/end point of alpha
		(bstart, bend): starting/end point of b
		nbf: number of frames
		path: you can specify a path where the files will be gererated (string ending with /)
		"""
		self.setup = setup
		self.V1s = V1s
		self.V1e = V1e
		self.Vzs = Vzs
		self.Vze = Vze
		self.nb = nb
		self.alpha_list = np.linspace(astart, aend, nbf)
		self.b_list = np.linspace(bstart, bend, nbf)
		self.path = path + "figs_for_movie/"

	def calc(self):
		
		os.popen("mkdir %s" % self.path)
		for i in range(len(self.alpha_list)):
			self.setup.alpha = self.alpha_list[i] 
			self.setup.b = self.b_list[i] 
			graph = Zgraph(self.setup,(self.V1s, self.V1e),(self.Vzs, self.Vze), self.nb)
			graph.calc()
			graph.plot()
			savefig("%simg%04d.png" % (self.path, i))
			del(graph)
	def mk_movie(self):
		os.popen("ffmpeg -qscale 1 -r 10 -b 9600 -i "+self.path+"img%04d.png movie.mp4")	
	def clean(self):
		os.popen("rm -R %s" % self.path)	
