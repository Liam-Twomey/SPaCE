import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.animation as animation
from scipy.spatial.transform import Rotation as R
from scipy.fft import fft,fftfreq,fftshift
import scipy.stats

class spin:
    def __init__(self,n:int,method:str='Uniform',width:float=1,dist_file:str=None,
        hyperfine:float=0, nuc:int=0,dipolar=0,S=1/2,ms=None):
        '''
        The ``spin`` class generates an object which represents a collection of spins

        Parameters
        ----------
        n : int
            Number of spins to generate 
        method: str
            'Uniform' to distribute the spins evenly across the bandwidth,
            'Histogram' to define a spectral shape (gaussian,etc.)
        width : float
            The bandwidth of the experiment, in GHz.
        dist_file: list
            list of lists, shape (2,n): [0] is the x axis of the spin population
            distribution, [1] is the  y-axis. Used only by the histogram mode.
        hyperfine: float 
            (endor only) value of the hyperfine coupling between the electron
            and the nucleus.
        nuc: 
            (endor only) number of nuclei coupled to the electron spin
        dipolar: int
            (endor only) the dipolar coupling between spins in the system (must be even)
        S: float
            (endor only) The electronic spin of the system.
        ms: str
            If not supplied, SPaCE assumes an even distribution across the m_s manifold.
            If a path to a Numpy binary file containing a 1D array of values between -S
            and S is provided, that population distribution is used instead.

        Attributes
        ----------
        n: int
            number of spins (set by ``n`` argument)
        nu: np.array

        ms: np.array
            Values of $m_s$ available to the spins.
        S: float
            electronic spin of the system
       

        hyperfine: float
            (ENDOR only) Hyperfine coupling value for these spins with an electron.
        n_nu: np.array
            (ENDOR only) Each spin's deviation from the Larmour frequency at the
            center of the band.
        n_ms: np.array
            (ENDOR only) Value of m_s for each spin

        pairing: np.array
            Used for distance measurements. (WIP)
        dipolar: np.array
            Dipolar coupling betweens spins. (WIP)
        
        Notes
        -----
        In the case of ENDOR mode, the n and nu attributes define the electron spin,
        while  n_nu and n_ms define the precession frequency and m_s manifold for the nuclei.
        '''
        self.n=n
        self.n_nu=np.zeros(n)
        self.hyperfine=hyperfine
        self.pairing=np.zeros(n)
        self.dipolar=np.zeros(n)
        self.S=S
        if ms==None:
            self.ms=np.arange(-S,S+.1,1)
        else:
            self.ms=np.load(ms)

        if method=='Uniform':
            self.nu=np.random.uniform(-width/2,width/2,n)
            self.n_ms=np.zeros(n)-1/2
        if method=='Histogram':
            if self.S==1/2:
                hist_data=np.histogram(dist_file[0],bins=1000,weights=dist_file[1])
                hist_dist = scipy.stats.rv_histogram(hist_data)
                self.nu=hist_dist.rvs(size=n)
                self.n_ms=np.zeros(n)-1/2
            else:
                prop=np.sum(dist_file[1:],axis=1)
                for i,j in enumerate(prop):
                    prop[i]=j/(self.S*(self.S+1)-self.ms[i]*(self.ms[i]+1))
                prop/=np.sum(prop)
                self.n_ms=np.random.choice(self.ms[:len(prop)],self.n,p=prop)
                self.nu=np.zeros(self.n)
                for i,j in enumerate(prop):
                    mask=self.n_ms==self.ms[i]
                    hist_data=np.histogram(dist_file[0],bins=1000,weights=dist_file[int(i+1)])
                    hist_dist=scipy.stats.rv_histogram(hist_data)
                    self.nu[mask]=hist_dist.rvs(size=np.sum(mask))

        if hyperfine!=0:
            dist=np.sin(np.linspace(0,np.pi,1000))
            hist_data=np.histogram(np.linspace(0,np.pi,1000),bins=1000,weights=dist)
            hist_dist=scipy.stats.rv_histogram(hist_data)
            angles=hist_dist.rvs(size=n)
            self.hyperfine=(3*np.cos(angles)**2-1)*hyperfine/2
            mask=angles>np.pi/2
            self.hyperfine[mask]*=-1
            self.n_nu=nuc+self.hyperfine

        if dipolar!=0:
            index=np.linspace(0,n-1,n)
            self.pairing=np.zeros_like(self.nu)
            if n % 2 !=0:
                return 'set n as an even number for dipolar coupling between spins'
            self.pairing[int(n/2):]=np.random.permutation(int(n/2))
            self.pairing[:int(n/2)]=np.argsort(self.pairing[int(n/2):])+int(n/2)
            dist=np.sin(np.linspace(0,np.pi,1000))
            hist_data=np.histogram(np.linspace(0,np.pi,1000),bins=1000,weights=dist)
            hist_dist=scipy.stats.rv_histogram(hist_data)
            angles=hist_dist.rvs(size=int(n/2))
            self.dipolar[:int(n/2)]=(3*np.cos(angles)**2-1)*dipolar/2
            self.dipolar[int(n/2):]=self.dipolar[:int(n/2)][np.argsort(self.pairing[:int(n/2)])]
    
    def __str__(self):
        return f"You have {self.n} spins"

class inversion:
    '''
    This class uses the fourier-transform method to calculate the inversion profile
    for a given pulse.

    Parameters
    ----------
    pulse: ``pulse`` object
        A RF pulse.

    '''

    def __init__(self,pulse):
        self.weights=fftshift(abs(fft(pulse.amp)))
        self.weights/=max(self.weights)
        self.freq=fftshift(fftfreq(pulse.time.shape[-1]))/pulse.step
        if pulse.type != 'rect':
            self.time=pulse.time
            self.f=self.time*((pulse.f1-pulse.f0)/pulse.tp)+pulse.f0

class trajectory:
    '''
    The precession sequence of the system

    Parameters
    ----------
    t: float
        Length of the precession period, in ns
    dt: float
        Time step to simulate at (temporal resolution)
    spin: ``spin`` object
        The spin system to precess.

    Attributes
    ----------
    
    ETC: misc
        All relevant parameters are imported from the ``spin`` object.

    Notes
    -----
    The ``trajectory`` will simulate t/dt time steps, separated by dt.
    '''
    def __init__(self,t,dt,spin):
        self.traj=np.zeros((spin.n,3,int(t//dt)))
        self.traj[:,2,0]=np.ones((spin.n))
        self.time=np.arange(0,t,dt)
        self.final_pos=self.traj[:,:,-1]
        self.S=spin.S
        self.n_ms=spin.n_ms
        self.nu=spin.nu
        self.n=spin.n
        self.dt=dt
        self.M=np.sum(self.traj,axis=0)
        self.deltaB=None
        self.angle=None
        self.n_nu=spin.n_nu
        self.hyperfine=spin.hyperfine
        self.pairing=spin.pairing
        self.tipped_angle=np.zeros(spin.n)
        self.dipolar=spin.dipolar
        self.contribution=np.ones(spin.n)
        self.traj_intensity=None

    def find_deltaB(self,central_nu=0):
        '''
        Calculates $/Delta B$ for the system.
        '''
        ge=2.0023
        be=9.27*10**-24 ##J/Tt
        h=(6.626*10**-34)*(10**9) ##J/Hz
        self.deltaB=(self.nu-central_nu)*h/ge/be
        
    def __str__(self):
        return f"trajectory for {self.time[-1]} ns"
    
    def precess(self,t0,t1):
        '''
        Calculates the precession between two timepoints.

        Parameters
        ----------
        t0 : float
            initial timepoint
        t1 : float
            final timepoint
        '''
        initial_pos=np.argmin(abs(self.time-t0))
        final_pos=np.argmin(abs(self.time-t1))
        for t in range(final_pos-initial_pos):
            phase=np.angle(self.traj[:,0,initial_pos+t]+self.traj[:,1,initial_pos+t]*1j) #find phase
            r=np.linalg.norm(self.traj[:,:2,initial_pos+t],axis=1)
            self.traj[:,0,initial_pos+t+1]=np.multiply(np.sin(2*np.pi*self.dt*self.nu+phase+np.pi/2),r) #update x
            self.traj[:,1,initial_pos+t+1]=np.multiply(np.cos(2*np.pi*self.dt*self.nu+phase+3*np.pi/2),r) #update y   
            self.traj[:,2,initial_pos+t+1]=self.traj[:,2,initial_pos+t]  
    
    def flip_angle(self,B1,tp):
        '''
        Calculates the flip angle.

        '''
        ge=2.0023
        be=9.27*10**-24 ##J/T
        h=6.626*10**-34/(2*np.pi) ##J/Hz/rad
        self.angle=ge*be*tp*10**-9/h*B1

    def flip(self,pulse,t,direction):
        '''
        Apply a ``pulse`` to this trajectory starting at time 't'.

        Parameters
        ----------
        pulse: ``pulse`` object
            A pulse to apply to the system
        t: float
            Time at which to apply the pulse
        direction: str
            The sign and axis of propagation of the pulse (i.e. +x, -y, etc.)
        '''
        if pulse.type=='rect':
            pos=np.argmin(abs(self.time-t))
            phases=['+x','-x','+y','-y']
            plane=[0,0,1,1]
            sign=[1,-1,1,-1]
            index=phases.index(direction)
            self.find_deltaB()
            theta=np.arctan(self.deltaB/pulse.B1)
            rotvec=np.zeros_like(self.traj[:,:,0])
            rotvec[:,plane[index]]=np.cos(theta)*sign[index]
            rotvec[:,2]=np.sin(theta)
            Beff=np.sqrt(np.multiply(self.deltaB,self.deltaB)+pulse.B1**2)
            self.flip_angle(Beff,self.dt)
            self.angle*=(self.S*(self.S+1)-np.multiply(self.n_ms,self.n_ms+1))**(1/2)
            for i in range(len(np.arange(t,t+pulse.tp,self.dt))):
                f_rotvec=np.multiply(self.angle,np.transpose(rotvec))
                r_flip=R.from_rotvec(np.transpose(f_rotvec))
                self.traj[:,:,pos+i+1]=r_flip.apply(self.traj[:,:,pos+i])
            if sum(self.dipolar)!=0:
                proportion=(self.traj[self.pairing,2,pos+i+1]+1)/2
                new_contribution=np.multiply(self.contribution,proportion)
                self.contribution=np.append(new_contribution,self.contribution-new_contribution)
                self.pairing=np.append(self.pairing,self.pairing)
                self.nu=np.append(self.nu,self.nu-self.dipolar)
                self.dipolar=np.append(self.dipolar,-1*self.dipolar[altered_spins.astype(int)])
                self.traj=np.append(self.traj,self.traj[altered_spins.astype(int),:,:],axis=0)
                self.contribution=np.append(self.contribution,-1*self.dipolar[altered_spins.astype(int)])

            # self.dipolar[altered_spins.astype(int)]=-1*self.dipolar[altered_spins.astype(int)]
            # self.nu[altered_spins.astype(int)]=self.nu[altered_spins.astype(int)]+self.dipolar[altered_spins.astype(int)]


        if pulse.type=='wurst':
            pos=np.argmin(abs(self.time-t))
            final_pos=np.argmin(abs(self.time-(t+pulse.tp)))
            p_traj=np.zeros((self.n,3,int(pulse.tp/pulse.resolution+1)))
            p_traj[:,:,0]=self.traj[:,:,pos]
            p_time=np.arange(t,t+pulse.tp,pulse.resolution)
            phases=['+x','-x','+y','-y']
            plane=[0,0,1,1]
            sign=[1,-1,1,-1]
            index=phases.index(direction)
            for i in range(len(np.arange(t,t+pulse.tp,pulse.resolution))):
                self.find_deltaB((i*pulse.resolution-pulse.tp/2)*pulse.BW/pulse.tp)
                CHIRP_amp=1-np.abs(np.sin(np.pi*(i*pulse.resolution-pulse.tp/2)/pulse.tp))**pulse.n
                theta=np.arctan(self.deltaB/(pulse.B1*CHIRP_amp))
                rotvec=np.zeros_like(p_traj[:,:,0])
                rotvec[:,plane[index]]=np.cos(theta)*sign[index]
                rotvec[:,2]=np.sin(theta)
                Beff=np.sqrt(np.multiply(self.deltaB,self.deltaB)+(pulse.B1*CHIRP_amp)**2)
                self.flip_angle(Beff,pulse.resolution)
                f_rotvec=np.multiply(self.angle,np.transpose(rotvec))
                r_flip=R.from_rotvec(np.transpose(f_rotvec))
                p_traj[:,:,i+1]=r_flip.apply(p_traj[:,:,i])
            index=[np.argmin(np.abs(p_time-x))+1 for x in np.arange(t,t+pulse.tp,self.dt)]
            self.traj[:,:,pos+1:final_pos+1]=p_traj[:,:,index]
            self.tipped_angle=np.arccos(self.traj[:,2,final_pos])
            self.nu=self.nu+np.multiply(self.dipolar,np.cos(self.tipped_angle))

        if pulse.type=='custom':
            pos=np.argmin(abs(self.time-t))
            final_pos=np.argmin(abs(self.time-(t+pulse.tp)))
            p_traj=np.zeros((self.n,3,int(pulse.tp/pulse.resolution+1)))
            p_traj[:,:,0]=self.traj[:,:,pos]
            p_time=np.arange(t,t+pulse.tp,pulse.resolution)
            phases=['+x','-x','+y','-y']
            plane=[0,0,1,1]
            sign=[1,-1,1,-1]
            index=phases.index(direction)
            for i in range(len(np.arange(t,t+pulse.tp,pulse.resolution))):
                self.find_deltaB(pulse.freq(i*pulse.resolution))
                amp=pulse.amp(i*pulse.resolution)
                theta=np.arctan(self.deltaB/(pulse.B1*amp))
                rotvec=np.zeros_like(p_traj[:,:,0])
                rotvec[:,plane[index]]=np.cos(theta)*sign[index]
                rotvec[:,2]=np.sin(theta)
                if pulse.B1*amp<0:
                    rotvec[:,plane[index]]*=-1
                    rotvec[:,2]*=-1
                Beff=np.sqrt(np.multiply(self.deltaB,self.deltaB)+(pulse.B1*amp)**2)
                self.flip_angle(Beff,pulse.resolution)
                f_rotvec=np.multiply(self.angle,np.transpose(rotvec))
                r_flip=R.from_rotvec(np.transpose(f_rotvec))
                p_traj[:,:,i+1]=r_flip.apply(p_traj[:,:,i])
            index=[np.argmin(np.abs(p_time-x))+1 for x in np.arange(t,t+pulse.tp,self.dt)]
            self.traj[:,:,pos+1:final_pos+1]=p_traj[:,:,index]
            self.tipped_angle=np.arccos(self.traj[:,2,final_pos])
            self.nu=self.nu+np.multiply(self.dipolar,np.cos(self.tipped_angle))

    def calc_intensity(self):
        '''
        Calculate the signal (echo) intensity of the trajectory.
        '''
        n_ms_mat=(self.S*(self.S+1)-self.n_ms*(self.n_ms+1))
        self.traj_intensity=np.einsum('i,ijk->ijk',n_ms_mat,self.traj)

    def net_M(self,ms=None):
        '''
        Calculate net [magnetization]?
        '''
        if ms==None:
            mask=self.n_ms!=0
        else:
            mask=self.n_ms==ms
        self.M=np.sum(self.traj_intensity[mask],axis=0)
        
    def flip_n(self,nu):
        '''
        Nuclear version of ``flip`` for use in ENDOR sim.
        '''
        mask=abs(self.n_nu-nu)<0.0003
        self.nu[mask]-=self.hyperfine[mask]

    def p_seq(self,x,t0,direction=None,norm_y=1,norm_x=1):
        '''
        Defines a pulse sequence in such a way it can be displayed easily,
        overlaid with the system's magnetization.

        Parameters
        ----------
        x: list
            list of ``pulse`` objects to add to the sequence
        t0: list of float
            list of starting timepoints for the pulses in ``x``
        direction: list of str
            list of directions along which to apply the pulses in ``x``, i.e.
            +x, -y, etc
        norm_y: float
            Factor by which to multiply the pulse intensity, to show it
            appropriately on the plot alongside magnetization.
        norm_x: float
            Factor by which to multiply the time axis of the pulse intensity,
            to show it appropriately on the plot overlaid with magnetization.
        '''
        self.seq=np.copy(self.M)
        for i,j in enumerate(x):
            index0=np.argmin(abs(self.time-t0[i]))
            index1=np.argmin(abs(self.time-(t0[i]+x[i].tp*norm_x)))
            if direction==None:
                self.seq[0,index0:index1]=np.ones_like(self.seq[0,index0:index1])*x[i].B1*norm_y
                self.seq[1,index0:index1]=np.zeros_like(self.seq[1,index0:index1])
            else:
                phases=['+x','-x','+y','-y']
                anti_phase=['+y','-y','+x','-x']
                plane=[0,0,1,1]
                sign=[1,-1,1,-1]
                self.seq[plane[phases.index(direction[i])],index0:index1]=np.ones_like(self.seq[plane[phases.index(direction[i])],
                                                                                                      index0:index1])*x[i].B1*norm_y*sign[phases.index(direction[i])]
                self.seq[plane[anti_phase.index(direction[i])],index0:index1]=np.zeros_like(self.seq[plane[anti_phase.index(direction[i])],
                                                                                                      index0:index1])
                
    def get_traj(self,nu):
        '''
        Return the trajectory.
        '''
        index=np.argmin(abs(self.nu-nu))
        return self.traj[index]
    
    def display_bloch(self,t0,t1,nu,filename,interval=400,writer='pillow'):
        '''
        Write out a gif of the bloch sphere representation of the trajectory.

        Parameters
        ----------
        t0: float
            Start time for animation
        t1 : float
            End time for animation
        nu: float
            Frequency (relative to Larmour) of spin to display
        filename: str
            Name of gif file to output
        interval: int
            Length which each frame is displayed (ms)
        writer: str
            A matplotlib.animation writer, dictates the engine used to print the output.
        '''

        fig = plt.figure()
        ax = plt.axes(projection='3d')

        if hasattr(nu,'__iter__'):
            s=[np.argmin(abs(self.nu-i)) for i in nu]
        else:
            s=np.argmin(abs(self.nu-nu))

        # Make sphere
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = np.outer(np.cos(u), np.sin(v))
        y = np.outer(np.sin(u), np.sin(v))
        z = np.outer(np.ones(np.size(u)), np.cos(v))

        n=int((t1-t0)/self.dt)
        ims = ['']*n
        start_pos=np.argmin(abs(self.time-t0))
        for i in range(n):
            if hasattr(s,'__iter__'):
                for j,k in enumerate(s):
                    if j==0:
                        xline=self.traj[k,0,start_pos:i+start_pos+1]
                        yline=self.traj[k,1,start_pos:i+start_pos+1]
                        zline=self.traj[k,2,start_pos:i+start_pos+1]
                        artist=ax.plot3D(xline, yline, zline, 'gray')
                    else:
                        xline=self.traj[k,0,start_pos:i+start_pos+1]
                        yline=self.traj[k,1,start_pos:i+start_pos+1]
                        zline=self.traj[k,2,start_pos:i+start_pos+1]
                        artist.append(ax.plot3D(xline, yline, zline, 'gray')[0])
                artist.append(ax.plot_surface(x, y, z,color='gray',alpha=0.1))
                t=np.linspace(0,50,51)
                xline=np.sin(t*np.pi*2/50)
                yline=np.cos(t*np.pi*2/50)
                zline=np.zeros_like(t)
                artist.append(ax.plot3D(xline, yline, zline, 'blue')[0])
                ims[i]=artist
            else:
                xline=self.traj[s,0,start_pos:i+start_pos+1]
                yline=self.traj[s,1,start_pos:i+start_pos+1]
                zline=self.traj[s,2,start_pos:i+start_pos+1]
                artist=ax.plot3D(xline, yline, zline, 'gray')
                artist.append(ax.plot_surface(x, y, z,color='gray',alpha=0.1))
                t=np.linspace(0,50,51)
                xline=np.sin(t*np.pi*2/50)
                yline=np.cos(t*np.pi*2/50)
                zline=np.zeros_like(t)
                artist.append(ax.plot3D(xline, yline, zline, 'blue')[0])
                ims[i]=artist
        # self.ims=ims
        ani = animation.ArtistAnimation(fig=fig, artists=ims, interval=interval)
        ani.save(filename=filename,writer=writer)

class pulse:
    '''
    An object representation of an RF/MW pulse.

    Parameters
    ----------
    None

    '''
    def __init__(self, alpha, tp, type='rect', f0=None, f1=None, n_wurst=None,
                 resolution=0.1,freq=None,amp=None,B1=None):
        self.type=type.lower()
        match self.type:
            case 'rect':
                assert((alpha and tp) is not  None)
                self.rect(alpha, tp)
            case 'wurst':
                assert((alpha and tp and f0 and f1 and n_wurst and resolution) is not  None)
                self.wurst(alpha,tp,f0,f1,n_wurst,resolution)
            case 'custom':
                assert((tp and freq and amp and B1 and resolution) is not None)
                self.custom(tp,freq,amp,B1,resolution)
            case _:
                return "Undefined type. Options are 'rect', 'wurst', or 'custom'."

    def __str__(self):
        match self.type:
            case 'rect': 
                return f"{self.tp} ns {self.type} pulse, with flip angle of {self.alpha/np.pi}*pi rad."
            case 'wurst': 
                return f"{self.tp} ns {self.type} pulse, with flip angle of {self.alpha/np.pi}*pi rad.\n \
                Pulse starts at {self.f0} GHz and ends at {self.f1} GHz, with a resolution of {self.resolution} GHz."
            case 'custom': 
                return f"{self.tp} ns pulse at {self.freq} GHz, with an amplitude of {self.amp}, with a B1 of {self.B1}." 

    def rect(self,alpha,tp):
        '''
        Generates a rectangular pulse of time tp and angle alpha. Used only by the class constructor.

        Parameters
        ----------
        alpha: float
            The pulse angle
        tp: float
            Length of the pulse (ns)
        '''
        # self.step=step
        self.tp=tp
        self.alpha = alpha
        ge=2.0023
        be=9.27*10**-24 ##J/T
        h=6.626*10**-34/(2*np.pi) ##J/Hz/rad
        self.B1=alpha/ge/be/(tp*10**-9)*h ##T
        # self.time=np.arange(0,time,step=step)
        # self.amp=np.zeros_like(self.time)
        # self.amp[np.argmin(abs(self.time-t0)):np.argmin(abs(self.time-(t0+tp)))]=np.ones_like(self.amp[np.argmin(abs(self.time-t0)):np.argmin(abs(self.time-(t0+tp)))])*amp
        #self.type='rect'
    
    def wurst(self,alpha,tp,f0,f1,n_wurst,resolution=0.1):
        '''
        Generates a WURST pulse of time tp and angle alpha, starting at f0 GHz and ending at f1 GHz. Used only by the class constructor.

        Parameters
        ----------
        alpha: float
            The pulse angle
        tp: float
            Length of the pulse (ns)
        f0: 
            Starting frequency
        f1: 
            Ending frequency
        n_wurst: 
            ??
        resolution:
            Frequency step??
        '''

        self.f0=f0
        self.f1=f1
        self.tp=tp
        self.alpha=alpha
        self.BW=np.abs(f1-f0)
        k=self.BW*10**9*2*np.pi/tp/10**-9 ## Hz**2
        ge=2.0023
        be=9.27*10**-24 ##J/T
        h=6.626*10**-34/(2*np.pi) ##J/Hz
        if alpha==np.pi:
            Q_crit=5
        if alpha==np.pi/2:
            Q_crit=2*np.log(2)/np.pi
        self.B1=np.sqrt(Q_crit*k)*h/ge/be
        self.n=n_wurst
        self.resolution=resolution
        #self.type='WURST'

    def custom(self,tp,freq,amp,B1,resolution=0.1):
        '''
        Method to define a custom pulse. Used only by the class constructor.

        Parameters
        ----------
        tp: float
            Pulse length
        freq: float
            Pulse frequency
        amp: float
            Pulse amplitude
        B1: float
            Applied magnetic field of pulse
        resolution: float
            Resolution for pulse calculation [units??]
        '''
        self.tp=tp
        self.freq=freq
        self.amp=amp
        self.B1=B1
        self.resolution=resolution
