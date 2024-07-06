import numpy as np
import h5py

class Dycore:
    def __init__(self, 
                 file: str):
        self.ds = h5py.File(file, "r")
        self.Rd = 287
        self.cp = 1004
        self.g  = 9.81
        self.H  = 6800
        # self.a  = 6.37122e6
        self.u  = self.getVar("grid_u_c_xyzt")
        # self.u_p  = self.getVar("grid_u_p_xyzt")
        # self.u_n  = self.getVar("grid_u_n_xyzt")
        
        self.v  = self.getVar("grid_v_c_xyzt")
        self.t  = self.getVar("grid_t_c_xyzt")
        self.ps = self.getVar("grid_ps_c_xyzt")
        self.p  = self.getVar("grid_p_full_xyzt")
        self.p_half  = self.getVar("grid_p_half_xyzt")

        # self.z_half  = self.getVar("grid_z_half_xyzt")
        self.z_full  = self.getVar("grid_z_full_xyzt")
        
        

        self.qv   = self.getVar("grid_tracers_c_xyzt")
        # self.qv_p = self.getVar("grid_tracers_p_xyzt")
        # self.qv_n = self.getVar("grid_tracers_n_xyzt")
        
        
        self.qv_diff = self.getVar("grid_tracers_diff_xyzt")
        self.w = self.getVar("grid_w_full_xyzt")
        
        # self.factor1 = self.getVar("factor1_xyzt")
        # self.factor2 = self.getVar("factor2_xyzt")
        # self.factor3 = self.getVar("factor3_xyzt")
        # self.factor4 = self.getVar("factor4_xyzt")


        # self.convection = self.getVar("convection_xyzt")
        
        
        ### the variables for plot
        self.sigma_mean           = np.nanmean(self.p/self.ps, axis=(0,3))
        self.sigma_onlyz          = np.nanmean(self.sigma_mean, axis=1)
        self.y                    = np.linspace(-90,90,64)
        self.yy, self.sigma_mean2 = np.meshgrid(self.y,self.sigma_onlyz)
        
        ### cooridate
        # self.x  = np.linspace(-180,180,128)
        # self.y  = np.linspace(-90,90,64)
        # self.xd = np.deg2rad(self.x)
        self.yd = np.deg2rad(self.y)
        # self.xx, self.yy = np.meshgrid(self.x,self.y)
        
        # WARNING: cos(-90) would be zero, doing this preventing it from divide by zero.
        self.cy     = np.cos(self.yd)
        self.cy[0]  = np.nan
        self.cy[-1] = np.nan

    def getVar(self, var):
        return np.asarray(self.ds[var])
    
    def cal_KE(self, u, v):
        self.KE = (u**2 + v**2)
        self.KE_mean = np.nanmean(self.KE, axis=(1,3))
        return self.KE_mean
    
    def cal_pre(self, qv_diff, p_half):
        g = 9.81
        Prec = np.zeros(qv_diff.shape)
        for i in range(1,20-1):
            Prec[:,i,:,:] = 1/g * qv_diff[:,i,:,:] * (p_half[:,i+1,:,:] - p_half[:,i,:,:])
        Prec[:, 0,:,:] = 1/g * qv_diff[:, 0,:,:] * (p_half[:, 1,:,:] - p_half[:, 0,:,:])
        # self.Prec[:,-1,:,:] = 1/self.g * self.qv[:,-1,:,:] * (self.p_half[:,-1,:,:] - self.p_half[:,-2,:,:])
        
        Prec_mean = np.nansum(Prec, axis=(1))
        Prec_mean2 = np.nanmean(Prec_mean, axis=(2))
        
        return Prec_mean2, Prec
    
    def cal_t(self):
        print(np.nanmax(self.t), np.nanmin(self.t))
        self.t_mean = np.nanmean(self.t[:,:,:,:], axis=(1,3))
        return self.t_mean
    
    def cal_t_last(self):
        print(np.nanmax(self.t), np.nanmin(self.t))
        self.t_mean = np.nanmean(self.t[-100:,:,:,:], axis=(1,3))
        return self.t_mean

