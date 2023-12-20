##<Usage: python Script.py {path_exp}>
### ex) python vaginal_pcr_update_reference.py "/home/kbkim/vaginal_pcr/input/EGvaginal_experiment_result.xlsx"

import os, datetime
import pandas as pd
import sys
import numpy as np
from scipy.stats import percentileofscore, pearsonr
import matplotlib.pyplot as plt

# Check if the script is being called with the correct arguments
if len(sys.argv) < 2:
    print("Usage: python Script.py <path_exp>")
    print("Example: python vaginal_pcr_update_reference.py \"/home/kbkim/vaginal_pcr/input/EGvaginal_experiment_result.xlsx\"")
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
class VaginalPCRUpdateRef:
    def __init__(self, path_exp, fplog=None):
        """
        Initializes a VaginalPCRUpdateRef object.

        Parameters:
        path_exp (str): Path of PCR experiment result file to analyze.
        """
        self.__fplog=fplog        
        
        ## Path of Reference files
        curdir = os.path.dirname(os.path.abspath(__file__))
        self.path_exp = path_exp
        self.path_db = f"{curdir}/input/EGvaginal_db_abundance.csv"
        
        ## Path of output files     


        
        ## Dataframe of Reference files
        self.df_exp = None
        self.df_db = None
        
        ## Dataframe of output files to calculate
        self.df_abundance = None
        
        ## Lists used for calculation
        self.li_new_sample_name = None
        self.li_microbiome = None
        
        
    def ReadDB(self):
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:           
            self.df_db = pd.read_csv(self.path_db)
            
            self.df_exp = pd.read_excel(self.path_exp)
            idx_well_id = list(self.df_exp[self.df_exp.columns[0]]).index('Well')
            self.df_exp = self.df_exp[idx_well_id:].rename(columns = self.df_exp[idx_well_id:].reset_index().iloc[0]).drop(idx_well_id).set_index('Well')
            
            li_column = self.df_exp.columns.values.tolist()
            
            if "Cт" in li_column:
                self.df_exp = self.df_exp[["Sample Name", "Target Name", "Cт"]]
                self.df_exp.rename(columns = {"Sample Name": "sample_name", "Target Name": "microbiome", "Cт": "Ct", "Tm1": "Tm1"}, inplace=True)
                
            elif "CT" in li_column:
                self.df_exp = self.df_exp[["Sample Name", "Target Name", "CT", 'Tm1']]
                self.df_exp.rename(columns = {"Sample Name": "sample_name", "Target Name": "microbiome", "CT": "Ct", "Tm1": "Tm1"}, inplace=True)
            
            else:
                print("Check the columns of Experiment result file")
                sys.exit()
            
            idx_nan_id = self.df_exp.index.to_list().index(np.nan)
            self.df_exp = self.df_exp.iloc[:idx_nan_id, :]
            self.df_exp.loc[self.df_exp['Ct'] == 'Undetermined', 'Ct'] = 40.1
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            
        return rv, rvmsg   

    def CalculateProportion(self):
        """
        Calculate the Relative Abundance and Save the Relative Abundance data as an Csv file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """           
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:      
            self.li_new_sample_name = list(dict.fromkeys(self.df_exp.sample_name)) 
            self.li_microbiome = ['L_crispatus', 'L_gasseri', 'L_iners', 'L_jensenii', 'G_vaginalis', 'F_vaginae', 'BVAB-1']
            
            self.df_abundance = pd.DataFrame(index = self.li_microbiome, columns = self.li_new_sample_name)
            self.df_abundance = self.df_abundance.fillna(0)
            
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_microbiome)):            
            
                    sample_name = self.li_new_sample_name[i]

                    condition_universal = (self.df_exp.sample_name == sample_name) & (self.df_exp.microbiome == 'Universal')
                    condition_microbiome = (self.df_exp.sample_name == sample_name) & (self.df_exp.microbiome == self.li_microbiome[j])

                    ct = self.df_exp[condition_microbiome]['Ct'].values[0]
                    ct_universal = self.df_exp[condition_universal]['Ct'].values[0]   
                    Tm1 = self.df_exp[condition_microbiome]['Tm1'].values[0]
                    
                    if (ct == 40.1) | (Tm1 <= 80):     
                        self.df_abundance.loc[self.li_microbiome[j], self.li_new_sample_name[i]] = 0
                        
                    else:  
                        self.df_abundance.loc[self.li_microbiome[j], self.li_new_sample_name[i]] = 2**(-(ct- ct_universal))
                        
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 
            
    # Insert data into DB - Merge the data frame df_db & df_exp
    def InsertDataDB(self): 
        """
        Inserts data into the database by merging the data frames df_db and df_abundance.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """   
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_db = pd.merge(self.df_db, self.df_abundance, how='outer',left_on='taxa', right_index=True, suffixes=['', '_right']) 
            self.df_db = self.df_db.fillna(0)
            self.df_db = self.df_db.filter(regex='^(?!.*_right).*') # Eliminate duplicate columns

            # Update the data - Convert df_exp to df_db
            self.df_abundance = self.df_db       
                      
            self.df_db_rev = self.df_db.set_index(keys=['taxa'], inplace=False, drop=True)    
            self.df_db_rev.to_csv(self.path_db, index_label='taxa')
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
            
        return rv, rvmsg      


####################################
# main
####################################
if __name__ == '__main__':
    
    vaginalupdate = VaginalPCRUpdateRef(path_exp)
    vaginalupdate.ReadDB()
    vaginalupdate.CalculateProportion()
    vaginalupdate.InsertDataDB()      
    
    print('Update Complete')     
            