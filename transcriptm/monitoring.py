#!/usr/bin/env python

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# CLASS: MONITORING
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class Monitoring:
    def __init__(self):
        self.reads = [0] * 6
        #{"raw","trimmed","phix","ncRNA","mapped","mapped_strict"}
   
   def get_tot_percentage(self):
        tot_percentage = [float(x)/self.reads[0]*100 for x in self.reads]
        return tot_percentage
    
    def get_percentage_prev(self):
        prev_percentage =[0] * 6
        prev_percentage[0]=100.0
        for i in range(1,6):
            prev_percentage[i] = round(float(self.reads[i]) /self.reads[i-1]*100,1)
        return prev_percentage

    