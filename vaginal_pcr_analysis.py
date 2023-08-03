##<Usage: python Script.py {path_exp}>
### ex) python vaginal_pcr_analysis.py "/home/kbkim/vaginal_pcr/input/experiment_result.xlsx"

import os, datetime
import pandas as pd
import sys
import numpy as np

# Check if the script is being called with the correct arguments
if len(sys.argv) < 2:
    print("Usage: python Script.py <path_exp>")
    print("Example: python vaginal_pcr_analysis.py \"/home/kbkim/vaginal_pcr/input/experiment_result.xls\"")
    sys.exit(1)
    
# path_exp : Path of Merged Proportion file to analyze
path_exp = sys.argv[1] 


#-------------------------------------------------------
# Common Function
#-------------------------------------------------------
def WriteLog(functionname, msg, type='INFO', fplog=None):
    #strmsg = "[%s][%s][%s] %s\n" % (datetime.datetime.now(), type, functionname, msg)
    #if( DEBUGMODE ): 
    #    print(strmsg)
    #else:
    #    if( fplog != None ):
    #        fplog.write(strmsg)
    
    head = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    writestr = f"[{head}][{functionname}] {msg}\n"
    #if( DEBUGMODE ):
    if( True ):
        #print(message)
        writestr = f"[{functionname}] {msg}\n"
        print(writestr)
        
    if( fplog != None ):
        fplog.write(writestr)
        fplog.flush()
        
###################################
# MainClass
###################################
class VaginalPCRAnalysis:
    def __init__(self, path_exp, fplog=None):
        """
        Initializes a VaginalPCRAnalysis object.

        Parameters:
        path_exp (str): Path of Merged Proportion file to analyze.
        """
        self.path_exp = path_exp
        self.__fplog=fplog        
        
        ## Path of Reference files
        curdir = os.path.dirname(os.path.abspath(__file__))

        
        ## Path of output files 
        self.path_harmful = f"{curdir}/output/harmful.csv"
        self.path_beneficial = f"{curdir}/output/beneficial.csv"        
        
        ## Dataframe of Reference files
        self.df_exp = None
        
        ## Dataframe of output files to calculate
        self.df_harmful = None
        self.df_beneficial = None        
        
        ## Lists used for calculation
        self.li_new_sample_name = None
        
        
        
    def ReadDB(self):
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:           
            self.df_exp = pd.read_excel(self.path_exp)
            idx_well_id = list(self.df_exp[self.df_exp.columns[0]]).index('Well')
            self.df_exp = self.df_exp[idx_well_id:].rename(columns = self.df_exp[idx_well_id:].reset_index().iloc[0]).drop(idx_well_id).set_index('Well')
            
            li_column = self.df_exp.columns.values.tolist()
            
            if "Cт" in li_column:
                self.df_exp = self.df_exp[["Sample Name", "Target Name", "Cт"]]
                self.df_exp.rename(columns = {"Sample Name": "sample_name", "Target Name": "microbiome", "Cт": "Ct"}, inplace=True)
                
            elif "CT" in li_column:
                self.df_exp = self.df_exp[["Sample Name", "Target Name", "CT"]]
                self.df_exp.rename(columns = {"Sample Name": "sample_name", "Target Name": "microbiome", "CT": "Ct"}, inplace=True)
            
            else:
                print("Check the columns of Experiment result file")
                sys.exit()
            
            idx_nan_id = self.df_exp.index.to_list().index(np.nan)
            self.df_exp = self.df_exp.iloc[:idx_nan_id, :]
            self.df_exp.loc[self.df_exp['Ct'] == 'Undetermined', 'Ct'] = 40.1
            print(self.df_exp)
                  
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            
        return rv, rvmsg   

    def CalculateProportion(self):
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:      
            self.li_new_sample_name = list(dict.fromkeys(self.df_exp.sample_name)) 
            
            self.df_exp['RA'] = 0
            
            for idx_exp, row_exp in self.df_exp.iterrows():
                
                sample_name = row_exp['sample_name']
                
                condition = (self.df_exp.sample_name == sample_name) & (self.df_exp.microbiome == 'Universal')
                
                ct = row_exp['Ct']
                ct_universal = self.df_exp[condition]['Ct'].values[0] 
                
                self.df_exp.loc[idx_exp, 'RA'] = 2**(-(ct- ct_universal))
            print(self.df_exp)
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            
####################################
# main
####################################
if __name__ == '__main__':
    
    vaginalpcranalysis = VaginalPCRAnalysis(path_exp)
    vaginalpcranalysis.ReadDB()      
    vaginalpcranalysis.CalculateProportion()  
    
    print('Analysis Complete')
