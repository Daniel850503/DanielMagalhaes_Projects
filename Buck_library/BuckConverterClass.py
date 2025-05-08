"""
Created on Wed May 8 2025

@author: daniel.magalhaes_xel
"""

class Converter_dcdc :
    def __init__(self, Vin_min, Vin_max, Vout_min, Vout_max, 
                 P_max, F_max, IL_ripple, Vout_ripple, Effmin):
            self.Vin_min = Vin_min
            self.Vin_max = Vin_max
            self.Vout_min = Vout_min
            self.Vout_max = Vout_max
            self.P_max = P_max
            self.F_max = F_max
            self.IL_ripple = IL_ripple # em percentagem
            self.Vout_ripple = Vout_ripple
            self.Effmin = Effmin
            
    def Ganho(self):
        self.Ganho_max = self.Vout_max/self.Vin_min
    
    def Out_current_dc (self):
        self.Iout_dc = round((self.P_max/self.Vout_min),2)
        
    def Out_current_max (self):    
        self.Iout_max = round((self.P_max/self.Vout_min)*
                               (1+((self.IL_ripple/100)/2)), 2)
     
    
    def Delta_current_max (self):
        self.DeltaIL_max = round((self.P_max/self.Vout_min)* 
                                 (self.IL_ripple/100), 2)
     
    def Out_current_min(self):
        self.Iout_min = round (self.Iout_max-self.DeltaIL_max,2)
        
        
    def total_Rdc_Loss_calc( self):
        self.Total_R_loss = round ((self.P_max * ( 1 - self.Effmin / 100 ))/
                                   (self.Iout_max*self.Iout_max))
    
    def total_Vdc_Loss_calc( self ):

        self.Total_V_loss = round ((self.P_max * ( 1 - self.Effmin / 100 ))/
                                   self.Iout_max, 2)
        
    def Loss_Max_calc (self):
        self.total_loss_estimation = round (self.P_max * 
                                    ( 1 - self.Effmin / 100 ),2)
    
    def duty_cycle_estimation_calc (self):
        self.duty_cycle_max = round((self.Vout_max + self.Total_V_loss) / 
                                    self.Vin_min, 2)
        
    def AC_DC_ratio_calc (self):
        self.current_ratio = round((self.DeltaIL_max / self.Iout_dc),2) 
    
    def Periode_calc (self):
        self.T_total = 1/self.F_max
    
    def Switch_time_on_calc (self):
        self.Ton = round (self.T_total * self.duty_cycle_max,9)
    
    def Switch_time_off_calc (self):
        self.Toff = round (self.T_total * (1-self.duty_cycle_max),9)
    
    def Current_ripple_Incrise_rate (self):
        self.Iinc_rate = self.DeltaIL_max / self.Ton
        
    def Current_ripple_decrise_rate (self):
        self.Idec_rate = self.DeltaIL_max / self.Toff