##<Usage: python Script.py {path_exp}>
### ex) python vaginal_pcr_analysis.py "/home/kbkim/vaginal_pcr/input/experiment_result.xlsx"

import os, datetime
import pandas as pd
import sys
import numpy as np
from scipy.stats import percentileofscore, pearsonr

# Check if the script is being called with the correct arguments
if len(sys.argv) < 2:
    print("Usage: python Script.py <path_exp>")
    print("Example: python vaginal_pcr_analysis.py \"/home/kbkim/vaginal_pcr/input/experiment_result.xlsx\"")
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
        path_exp (str): Path of PCR experiment result file to analyze.
        """
        self.__fplog=fplog        
        
        ## Path of Reference files
        curdir = os.path.dirname(os.path.abspath(__file__))
        self.path_exp = path_exp
        self.path_db = f"{curdir}/input/db_abundance.csv"
        
        ## Path of output files     
        self.path_percentile_rank_output = f"{curdir}/output/percentile_rank.csv"
        self.path_eval_output = f"{curdir}/output/eval.csv"
        
        ## Dataframe of Reference files
        self.df_exp = None
        self.df_db = None
        
        ## Dataframe of output files to calculate
        self.df_abundance = None
        self.df_percentile_rank = None 
        self.df_eval = None
        
        ## Lists used for calculation
        self.li_new_sample_name = None
        self.li_microbiome = None
        
        
    def ReadDB(self):
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:           
            self.df_db = pd.read_csv(self.path_db, index_col=0)
            
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
                    
                    self.df_abundance.loc[self.li_microbiome[j], self.li_new_sample_name[i]] = 2**(-(ct- ct_universal))
            #self.df_abundance = self.df_abundance.rename_axis('taxa', axis=1)                        
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process") 

    def CalculatePercentileRank(self):
        """
        Calculate the Percentile Rank and Save the Percentile Rank data as an Csv file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:      

            # Create an empty data frame 
            self.df_percentile_rank = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_microbiome)
                              
            # Loop through all samples and phenotypes and calculate the percentile rank
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_microbiome)):
                    self.df_percentile_rank.loc[self.li_new_sample_name[i], self.li_microbiome[j]] = (percentileofscore(list(self.df_db.loc[self.li_microbiome[j]]), self.df_abundance.loc[self.li_microbiome[j], self.li_new_sample_name[i]], kind='mean')).round(1)
            

            # Outliers
            # Replace percentile ranks that are less than or equal to 5 with 5, and those that are greater than or equal to 95 with 95
            for i in range(len(self.li_microbiome)):
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_microbiome[i]]<=5, self.li_microbiome[i]] = 5.0
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_microbiome[i]]>=95, self.li_microbiome[i]] = 95.0      
                     
            # Replace missing values with the string 'None'    
            self.df_percentile_rank = self.df_percentile_rank.fillna('None')

            # Save the output file - Percentile Rank of the samples
            self.df_percentile_rank.to_csv(self.path_percentile_rank_output, encoding="utf-8-sig", index_label='serial_number')
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg    
    
    def EvaluatePercentileRank(self):
        """
        Evaluate based on percentile rank value

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:                 
            self.df_eval = pd.DataFrame(index=self.df_percentile_rank.index)
            
            li_positive_var = ['L_crispatus', 'L_gasseri', 'L_iners', 'L_jensenii']
            
            for col in self.df_percentile_rank:
                if col in li_positive_var:
                    # Define the conditions and corresponding values
                    conditions = [
                        self.df_percentile_rank[col] >= 70,
                        (self.df_percentile_rank[col] >= 30) & (self.df_percentile_rank[col] < 70),
                        self.df_percentile_rank[col] <= 30
                    ]
                    
                    values = ['좋음', '보통', '나쁨']     
                    
                    self.df_eval[col] = np.select(conditions, values)  
                                 
                else: 
                    # Define the conditions and corresponding values
                    conditions = [
                        self.df_percentile_rank[col] >= 70,
                        (self.df_percentile_rank[col] >= 30) & (self.df_percentile_rank[col] < 70),
                        self.df_percentile_rank[col] <= 30
                    ]
                    
                    values = ['나쁨', '보통', '좋음']     
                    
                    self.df_eval[col] = np.select(conditions, values)       

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg         

    def ClassificateType(self):
        """
        Classificate a Type based on percentile rank value and Save the Evaluation data as an Csv file

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:  
            dict_score = self.df_percentile_rank.to_dict('records')
            
            dict_type = {'L_crispatus': 'CST I', 'L_gasseri': 'CST II', 'L_iners': 'CST III',  
                         'G_vaginalis': 'CST IV-A', 'F_vaginae': 'CST IV-B', 'BVAB-1': 'CST IV-C', 'L_jensenii': 'CST V'}

            for i in range(len(self.li_new_sample_name)):     
                max_taxa = max(dict_score[i],key=dict_score[i].get)
                
                self.df_eval.loc[self.li_new_sample_name[i], 'Type'] = dict_type[max_taxa]
                    
            # Save the output file - df_eval
            self.df_eval.to_csv(self.path_eval_output, encoding="utf-8-sig", index_label='serial_number') 

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
    
    vaginalpcranalysis = VaginalPCRAnalysis(path_exp)
    vaginalpcranalysis.ReadDB()      
    vaginalpcranalysis.CalculateProportion()  
    vaginalpcranalysis.CalculatePercentileRank()   
    vaginalpcranalysis.EvaluatePercentileRank()     
    vaginalpcranalysis.ClassificateType()    
    
    print('Analysis Complete')
