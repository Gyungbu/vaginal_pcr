import os, datetime
import pandas as pd
import sys
import numpy as np
from scipy.stats import percentileofscore, pearsonr
import matplotlib.pyplot as plt

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
        
# Histogram Plot 
def save_histograms_to_file(df, filename):
    num_rows = df.shape[1]
    fig, axs = plt.subplots(num_rows, 1, figsize=(8, 6*num_rows))
    
    for i in range(num_rows):           
        data = df.iloc[:,i]
        axs[i].hist(data, weights=np.ones(len(data)) / len(data)*100, bins=[0,20,40,60,80,100])
        axs[i].set_title(df.columns.to_list()[i][:-9])
        axs[i].set_xlabel('Sum of Relative Abundance[%]')
        axs[i].set_ylabel('Percentage of samples[%]')        
        axs[i].set_xlim([0, 100])
        
        if i == 0:
            beneficial_distribution = axs[i].hist(data, weights=np.ones(len(data)) / len(data)*100, bins=[0,20,40,60,80,100])[0]
        if i == 1:
            harmful_distribution = axs[i].hist(data, weights=np.ones(len(data)) / len(data)*100, bins=[0,20,40,60,80,100])[0]        
    
    str_beneficial_distribution = ','.join(str(x) for x in beneficial_distribution)
    str_harmful_distribution = ','.join(str(x) for x in harmful_distribution)

    plt.tight_layout()
    plt.savefig(filename)          
    
    return str_beneficial_distribution, str_harmful_distribution
    
    
###################################
# MainClass
###################################
class VaginalPCRAnalysis:
    def __init__(self, path_exp, outdir=None, fplog=None):
        """
        Initializes a VaginalPCRAnalysis object.

        Parameters:
        path_exp (str): Path of PCR experiment result file to analyze.
        """
        self.__fplog=fplog        
        
        ## Path of Reference files
        curdir = os.path.dirname(os.path.abspath(__file__))
        self.path_exp = path_exp
        self.path_db = f"{curdir}/input/EGvaginal_db_abundance.csv"
                       
        ###output
        if( outdir is not None ):
            self.outdir = outdir
        else:
            self.outdir = f"{curdir}/output"         
            
            
        ## Path of output files     
        self.path_abundance_output = f"{self.outdir}/EGvaginal_abundance.csv"
        self.path_eval_output = f"{self.outdir}/EGvaginal_eval.csv"
        self.path_mean_abundance = f"{self.outdir}/EGvaginal_mean_abundance.csv"
        self.path_hist = f"{self.outdir}/EGvaginal_abundance_hist.png"
        
        ## Dataframe of Reference files
        self.df_exp = None
        self.df_db = None
        
        ## Dataframe of output files to calculate
        self.df_abundance = None
        self.df_eval = None
        self.df_mean_abundance = None 
        
        ## Lists used for calculation
        self.li_new_sample_name = None
        self.li_microbiome = None
        
        ## Dictionaries used for calculation        
        self.dict_mean_abundance = None
        
        
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
                self.df_exp = self.df_exp[["Sample Name", "Target Name", "Cт", 'Tm1']]
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
            self.df_abundance = self.df_abundance.rename_axis('serial_number', axis=1)                        
            
            self.df_abundance = self.df_abundance.transpose()
            
            # Save the output file - Abundance of the samples
            self.df_abundance.to_csv(self.path_abundance_output, encoding="utf-8-sig", index_label='serial_number')      
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")   
    
    def EvaluateProportion(self):
        """
        Evaluate based on proportion value

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:                 
            self.df_eval = pd.DataFrame(index=self.df_abundance.index)

            self.df_mean_abundance = self.df_db.mean(axis=1, numeric_only=True).to_frame()     
            self.df_mean_abundance.columns =['value']
            
            self.dict_mean_abundance = self.df_db.mean(axis=1, numeric_only=True).to_dict()
            
            
            
            for col in self.df_abundance:
                # Define the conditions and corresponding values
                conditions = [
                    self.df_abundance[col] >= 0.5,
                    (self.df_abundance[col] >= self.dict_mean_abundance[col]) & (self.df_abundance[col] < 0.5),
                    self.df_abundance[col] < self.dict_mean_abundance[col]
                ]

                values = ['높음', '보통', '낮음']     

                self.df_eval[col] = np.select(conditions, values)  
                                 

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg         

    def ClassifyType(self):
        """
        Classify a Type based on percentile rank value and Save the Evaluation data as an Csv file

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:  
            dict_abundance = self.df_abundance.to_dict('records')
            dict_type = {'L_crispatus': '항균든든', 'L_gasseri': '항균유지', 'L_iners': '면역주의', 'L_jensenii': '항균특별',  
                         'G_vaginalis': '면역저하', 'F_vaginae': '면역저하', 'BVAB-1': '면역저하'}
            
            self.df_eval['SprayType'] = 'C'
            
            for idx in range(len(self.li_new_sample_name)):                  
                total_abundance = sum(dict_abundance[idx].values())
                
                sum_beneficial = sum(list(dict_abundance[idx].values())[0:4])
                sum_harmful = sum(list(dict_abundance[idx].values())[4:])
                
                if total_abundance < 0.05:
                    self.df_eval.loc[self.li_new_sample_name[idx], 'Type'] = '기타유형'
                    self.df_eval.loc[self.li_new_sample_name[idx], 'SprayType'] = 'G'
                    
                elif sum_harmful > sum_beneficial:
                    self.df_eval.loc[self.li_new_sample_name[idx], 'Type'] = '면역저하'
                    self.df_eval.loc[self.li_new_sample_name[idx], 'SprayType'] = 'G'                    
                    
                else:   
                    max_taxa = max(dict_abundance[idx],key=dict_abundance[idx].get)
                    self.df_eval.loc[self.li_new_sample_name[idx], 'Type'] = dict_type[max_taxa]
                    
                    if dict_type[max_taxa] == '면역저하':
                        self.df_eval.loc[self.li_new_sample_name[idx], 'SprayType'] = 'G'                    
                    
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg         
    
    def CalculateTotalAbundance(self):
        """
        Classify a Type based on percentile rank value and Save the Evaluation data as an Csv file

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:  
            self.df_eval['beneficial_total[%]'] = (self.df_abundance['L_crispatus'] + self.df_abundance['L_gasseri'] + self.df_abundance['L_iners'] + self.df_abundance['L_jensenii'])*100
              
            self.df_eval['harmful_total[%]'] = (self.df_abundance['G_vaginalis'] + self.df_abundance['F_vaginae'] + self.df_abundance['BVAB-1'])*100 
            
            self.df_db = self.df_db.transpose() 
            self.df_db['beneficial_total[%]'] = (self.df_db['L_crispatus'] + self.df_db['L_gasseri'] + self.df_db['L_iners'] + self.df_db['L_jensenii'])*100
              
            self.df_db['harmful_total[%]'] = (self.df_db['G_vaginalis'] + self.df_db['F_vaginae'] + self.df_db['BVAB-1'])*100                         
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg     

    def EvaluateBeneficialHarmful(self):
        """
        Evaluate Beneficial & Harmful microbiome

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:                                                    
            self.dict_mean_abundance = self.df_db.mean(axis=0, numeric_only=True).to_dict()          
            
            for col in ['beneficial_total[%]', 'harmful_total[%]']:
                # Define the conditions and corresponding values
                conditions = [
                    self.df_eval[col] >= 50,
                    (self.df_eval[col] >= self.dict_mean_abundance[col]) & (self.df_eval[col] < 50),
                    self.df_eval[col] < self.dict_mean_abundance[col]
                ]

                values = ['높음', '보통', '낮음']     

                self.df_eval[f"{col[:-3]}_eval"] = np.select(conditions, values)  
            # Save the output file - df_eval
            self.df_eval.to_csv(self.path_eval_output, encoding="utf-8-sig", index_label='serial_number')                          
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg          
    
    def PlotDistribution(self): 
        """
        Plot the Distribution - Relative Abundance of Harmful & Beneficial microbiome  

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """   
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        rv = True
        rvmsg = "Success"
        
        try:          
            # Histogram Plot - mrs 
            save_histograms_to_file(self.df_db[['beneficial_total[%]', 'harmful_total[%]']], self.path_hist)
                
            self.df_mean_abundance.loc['beneficial_distribution'] = save_histograms_to_file(self.df_db[['beneficial_total[%]', 'harmful_total[%]']], self.path_hist)[0]
            self.df_mean_abundance.loc['harmful_distribution'] = save_histograms_to_file(self.df_db[['beneficial_total[%]', 'harmful_total[%]']], self.path_hist)[1]
                  
            # Save the output file - Abundance of the samples
            self.df_mean_abundance.to_csv(self.path_mean_abundance, encoding="utf-8-sig", index_label='taxa')       
            
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
    
    path_exp = "input/EGvaginal_experiment_result.xlsx"
    
    vaginalpcranalysis = VaginalPCRAnalysis(path_exp)
    vaginalpcranalysis.ReadDB()      
    vaginalpcranalysis.CalculateProportion()   
    vaginalpcranalysis.EvaluateProportion()     
    vaginalpcranalysis.ClassifyType() 
    vaginalpcranalysis.CalculateTotalAbundance()     
    vaginalpcranalysis.EvaluateBeneficialHarmful() 
    vaginalpcranalysis.PlotDistribution()    

    
    print('Analysis Complete')
