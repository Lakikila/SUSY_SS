#define Nova_cxx
#include "Nova.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <algorithm>

void SetActiveBranches(TTree *fChain){ //setting up the branches that we need
	fChain->SetBranchStatus("*",0);
   	fChain->SetBranchStatus("truth_NeutralinoFromGluino_eta",1);
   	fChain->SetBranchStatus("truth_NeutralinoFromGluino_phi",1);
   	fChain->SetBranchStatus("truth_NeutralinoFromGluino_pt",1);
   	fChain->SetBranchStatus("truth_NeutralinoFromGluino_e",1);

   	fChain->SetBranchStatus("truth_QuarkFromNeutralino_eta",1);
   	fChain->SetBranchStatus("truth_QuarkFromNeutralino_phi",1);
   	fChain->SetBranchStatus("truth_QuarkFromNeutralino_pt",1);
   	fChain->SetBranchStatus("truth_QuarkFromNeutralino_e",1);

	fChain->SetBranchStatus("truth_QuarkFromGluino_eta",1);
   	fChain->SetBranchStatus("truth_QuarkFromGluino_phi",1);
   	fChain->SetBranchStatus("truth_QuarkFromGluino_pt",1);
   	fChain->SetBranchStatus("truth_QuarkFromGluino_e",1);

   	fChain->SetBranchStatus("fatjet_eta",1);
   	fChain->SetBranchStatus("fatjet_phi",1);
   	fChain->SetBranchStatus("fatjet_pt",1);
   	fChain->SetBranchStatus("fatjet_e",1);

	fChain->SetBranchStatus("jet_eta",1);
   	fChain->SetBranchStatus("jet_phi",1);
   	fChain->SetBranchStatus("jet_pt",1);
   	fChain->SetBranchStatus("jet_e",1);

	fChain->SetBranchStatus("fatjet_tau1_wta",1);
	fChain->SetBranchStatus("fatjet_tau21_wta",1);
	fChain->SetBranchStatus("fatjet_tau2_wta",1);
	fChain->SetBranchStatus("fatjet_tau32_wta",1);
	fChain->SetBranchStatus("fatjet_tau3_wta",1);
	fChain->SetBranchStatus("fatjet_ungrtrk500",1);
	fChain->SetBranchStatus("fatjet_NTrimSubjets",1);
	

}
void printTLorentz(TLorentzVector v){ //tlorentz printer
	cout<< "TLorentz: (" << v.Pt() << ", "<< v.Eta()<< ", " << v.Phi()<< ", " << v.E()<< ")"<<endl;
}
void settingNeutralinos(vector<TLorentzVector>& vn,vector<double> *pt,vector<double> *eta,vector<double> *phi,vector<double> *e){
	//setting neutralinos up at beginning of every process   	
	vn=vector<TLorentzVector>();
	TLorentzVector v;
	v.SetPtEtaPhiE((*pt)[0],(*eta)[0],(*phi)[0],(*e)[0]);
	vn.push_back(v);
	v.SetPtEtaPhiE((Double_t)(*pt)[1],(Double_t)(*eta)[1],(Double_t)(*phi)[1],(Double_t)(*e)[1]);
	vn.push_back(v);
	//there's no need for (Double_t) cast
}
void findingMin(vector<TLorentzVector>& vn,vector<double> *pt,vector<double> *eta,vector<double> *phi,vector<double> *e,vector<TLorentzVector>& min_fatjet,vector<TLorentzVector>& fatjet_rest){ 	//function that puts the fatjat that has minimal dR with his neutralino.
		//vector_neutralino[i] and his fatjet is min_fatjet[i]
      	TLorentzVector fatjet;
      	fatjet.SetPtEtaPhiE((*pt)[0],(*eta)[0],(*phi)[0],(*e)[0]);
      	double minimum=vn[0].DeltaR(fatjet);

      	int index_min=0;
      	int index_n=0;
      	for(int i=0; i< pt->size();i++){
          	fatjet.SetPtEtaPhiE((*pt)[i],(*eta)[i],(*phi)[i],(*e)[i]);
          	if(vn[0].DeltaR(fatjet) < vn[1].DeltaR(fatjet)){
              		if( vn[0].DeltaR(fatjet)<=minimum){
				minimum=vn[0].DeltaR(fatjet);
                  		index_min=i;
                  		index_n=0;
                  		min_fatjet[0]=fatjet;
              		}
          	}else{
              		if(vn[1].DeltaR(fatjet)<=minimum){
                  		minimum=vn[1].DeltaR(fatjet);
                  		index_min=i;
                  		index_n=1;
                  		min_fatjet[1]=fatjet;
              		}	
         	} 
      } //first jet is set up and his neutralino is vector_neutralino[index_n]
	//index of first jet among other fatjets is index_min
      index_n=(index_n+1)%2;
      int index_help=(index_min+1)%(pt->size()-1);
      fatjet.SetPtEtaPhiE((*pt)[index_help],(*eta)[index_help],(*phi)[index_help],(*e)[index_help]);
      
      minimum=vn[index_n].DeltaR(fatjet);
      for(int i=0; i< pt->size();i++){
          	if(i!=index_min){ 	//we ignore the jet that is mathced to previous neutralino
		    	fatjet.SetPtEtaPhiE((*pt)[i],(*eta)[i],(*phi)[i],(*e)[i]);
            		if(vn[index_n].DeltaR(fatjet) <= minimum ){
                	  	minimum=vn[index_n].DeltaR(fatjet);
                	  	min_fatjet[index_n]=fatjet;
            		}
      		}
      }
      for(int i=0; i< pt->size(); i++){
		if( (*pt)[i] != min_fatjet[0].Pt() && (*pt)[i] != min_fatjet[1].Pt() ){
            	fatjet.SetPtEtaPhiE((*pt)[i],(*eta)[i],(*phi)[i],(*e)[i]);
		fatjet_rest.push_back(fatjet);
		}
      }
}

bool quark_from_gluino_matching(vector<TLorentzVector>& quarks,vector<double> *pt,vector<double> *eta,vector<double> *phi,vector<double> *e,vector<TLorentzVector>& jet_pass){
	
	vector<TLorentzVector> copy_jet;
	vector<TLorentzVector> copy_quarks;
	vector<double> minimumi;
	double tmp_minimum=500,dR_tmp;
	
	bool break_flag=false;

	

	quarks.resize(pt->size());

	for(int i=0;i< pt->size();i++){
		quarks[i].SetPtEtaPhiE((*pt)[i],(*eta)[i],(*phi)[i],(*e)[i]);
	}

	copy_quarks=quarks;
	copy_jet= jet_pass;
	if(copy_jet.size()<4){
		return false;
	}
	cout<< quarks.size()<< " " <<copy_quarks.size()<<endl;
	cout<<jet_pass.size()<<" "<<copy_jet.size()<<endl;
	
	for(int i=0;i<copy_quarks.size();i++){
		for(int j=0;j<copy_jet.size();j++){
			dR_tmp= (copy_quarks)[i].DeltaR((copy_jet)[j]);
			cout << dR_tmp <<endl;
			if(dR_tmp < tmp_minimum ){
				tmp_minimum=dR_tmp;
				
			} 
		}
	}
	
	for(auto it1= copy_quarks.begin(); it1!= copy_quarks.end(); it1++ ){
		for(auto it2= copy_jet.begin(); it2!= copy_jet.end(); it2++ ){
			dR_tmp=(*it1).DeltaR(*it2);
			if(dR_tmp == tmp_minimum ){
				
				
				copy_jet.erase(it2);
				copy_quarks.erase(it1);
				
				break_flag=true;
				break;
				
			}
			
		}
		if(break_flag==true){
				break_flag=false;
				break;
		}
	}
	minimumi.push_back(tmp_minimum);
	
	cout <<  " a njihova razlika je " << tmp_minimum << endl;
	tmp_minimum=500;
	cout<< quarks.size()<< " " <<copy_quarks.size()<<endl;
	cout<<jet_pass.size()<<" "<<copy_jet.size()<<endl;
	for(int i=0;i<copy_quarks.size();i++){
		for(int j=0;j<copy_jet.size();j++){
			dR_tmp=copy_quarks[i].DeltaR(copy_jet[j]);
			if(dR_tmp < tmp_minimum ){
				tmp_minimum=dR_tmp;
				
			}
		}
	}
	
	for(auto it1= copy_quarks.begin(); it1!= copy_quarks.end(); it1++ ){
		for(auto it2= copy_jet.begin(); it2!= copy_jet.end(); it2++ ){
			dR_tmp=(*it1).DeltaR(*it2);
			if(dR_tmp == tmp_minimum ){
				
				
				copy_jet.erase(it2);
				copy_quarks.erase(it1);
				
				break_flag=true;
				break;
				
			}
			
		}
		if(break_flag==true){
				break_flag=false;
				break;
		}
	}
	minimumi.push_back(tmp_minimum);
	
	cout << " a njihova razlika je " << tmp_minimum << endl;
	tmp_minimum=500;
	cout<< quarks.size()<< " " <<copy_quarks.size()<<endl;
	cout<<jet_pass.size()<<" "<<copy_jet.size()<<endl;	
	for(int i=0;i<copy_quarks.size();i++){
		for(int j=0;j<copy_jet.size();j++){
			dR_tmp=copy_quarks[i].DeltaR(copy_jet[j]);
			if(dR_tmp < tmp_minimum ){
				tmp_minimum=dR_tmp;
				
			}
		}
	}
	
	for(auto it1= copy_quarks.begin(); it1!= copy_quarks.end(); it1++ ){
		for(auto it2= copy_jet.begin(); it2!= copy_jet.end(); it2++ ){
			dR_tmp=(*it1).DeltaR(*it2);
			if(dR_tmp == tmp_minimum ){
				
				
				
				copy_jet.erase(it2);
				copy_quarks.erase(it1);
				
				break_flag=true;
				break;
				
			}
			
		}
		if(break_flag==true){
				break_flag=false;
				break;
		}
	}
	minimumi.push_back(tmp_minimum);
	
	cout << " a njihova razlika je " << tmp_minimum << endl;
	tmp_minimum=500;
	cout<< quarks.size()<< " " <<copy_quarks.size()<<endl;
	cout<<jet_pass.size()<<" "<<copy_jet.size()<<endl;	
	for(int i=0;i<copy_quarks.size();i++){
		for(int j=0;j<copy_jet.size();j++){
			dR_tmp=copy_quarks[i].DeltaR(copy_jet[j]);
			if(dR_tmp < tmp_minimum ){
				tmp_minimum=dR_tmp;
				
			}
		}
	}
	
	for(auto it1= copy_quarks.begin(); it1!= copy_quarks.end(); it1++ ){
		for(auto it2= copy_jet.begin(); it2!= copy_jet.end(); it2++ ){
			dR_tmp=(*it1).DeltaR(*it2);
			if(dR_tmp == tmp_minimum ){
				
				
				
				copy_jet.erase(it2);
				copy_quarks.erase(it1);
				
				break_flag=true;
				break;
				
			}
			
		}
		if(break_flag==true){
				break_flag=false;
				break;
		}
	}
	minimumi.push_back(tmp_minimum);
	cout << " a njihova razlika je " << tmp_minimum << endl;
	tmp_minimum=500;
	for(int i=0;i< minimumi.size();i++){
		if(minimumi[i] > 0.8 ){
			return false;
		}
	}
	
	return true;	
	
}
			

void CheckIfQuarkIsOutOfJet(TLorentzVector quark, TLorentzVector jet,int& tmp_lost_quarks,int& total_lost_quarks,bool& process_unsuccessful_flag){
	//test if quark and jet are mismatched	
	if(quark.DeltaR(jet) > 0.8){
		process_unsuccessful_flag=true;
		tmp_lost_quarks++;
		total_lost_quarks++;
	}
}
void checkingForReasonOfProcessFailure(int& q,int& nn,int& nj,int& q_nn,int& q_nj,int& nn_nj,int& q_nn_nj,int& u,
		int& tmp_lost_quarks,vector<TLorentzVector> vn,vector<TLorentzVector>  min_fatjet){
	//investigation why the process is unsuccessful
	bool q_flag= (tmp_lost_quarks==1)?true:false;		//q_flag indicates that only 1 quark failed in process
	bool nn_flag= (vn[0].DeltaR(vn[1]) <= 1)?true:false;	//nn_flag indicates that neutralinos in process are too close to each other
	bool nj_flag=(vn[0].DeltaR(min_fatjet[0])>=0.8 || vn[1].DeltaR(min_fatjet[1])>=0.8)?true:false; //nj_flag indicates that we lost neutralino jet
	if(!q_flag && !nn_flag && !nj_flag) u++;
	if(!q_flag && !nn_flag && nj_flag) nj++;
	if(!q_flag && nn_flag && !nj_flag) nn++;
	if(q_flag && !nn_flag && !nj_flag) q++;
	if(!q_flag && nn_flag && nj_flag) nn_nj++;
	if(q_flag && !nn_flag && nj_flag) q_nj++;
	if(q_flag && nn_flag && !nj_flag) q_nn++;
	if(q_flag && nn_flag && nj_flag) q_nn_nj++;

}
void resetFlags(bool& process_unsuccessful_flag, int& tmp_lost_quarks, vector<TLorentzVector>& min_fatjet){
	//resetting variables for next process	
	process_unsuccessful_flag=false;	
	tmp_lost_quarks=0;
	min_fatjet=vector<TLorentzVector>(2);
}	
void myPrint(int total_lost_quarks,int process_unsuccessful,Long64_t nentries,int nn,int nj,int q,int nn_nj,int q_nn,int q_nj,int q_nn_nj,int u){
	//printing our results
    	cout<< "# of unmatched quarks: " << total_lost_quarks << endl;
    	cout<< "# of unsuccessful processes: "<< process_unsuccessful<< " and its % is: " << process_unsuccessful*1.0/nentries*100<<"%" << endl;
    	cout<< "out of all unsuccessful proccesses:\n";
	cout << "\t-(pure) neutralinos are too close in "<< nn << " and that is " << nn*1.0/process_unsuccessful*100 <<"%" << endl;
    	cout<< "\t-(pure) at least one neutralino fatjet is lost in "<< nj << " and that is " << nj*1.0/process_unsuccessful*100 <<"%" << endl;
    	cout<< "\t-(pure) a single quark mismatched in "<< q << " and that is " << q*1.0/process_unsuccessful*100 <<"%" << endl;
    	cout<< "\t-neutralinos too close & at least one neutralino fatjet lost in "<< nn_nj<< 
		" and that is " << nn_nj*1.0/process_unsuccessful*100 <<"%" << endl;
    	cout<< "\t-neutralinos are too close & exactly one quark is mismatched in "<< q_nn << 
		" and that is " << q_nn*1.0/process_unsuccessful*100 <<"%" << endl;
    	cout<< "\t-at least one neutralino jet is lost & exactly one quark is mismatched "<< q_nj << 
		" and that is " << q_nj*1.0/process_unsuccessful*100 <<"%" << endl;
    	cout<< "\t-neutralinos are too close & at least one neutralino fatjet lost & exactly one quark is mismatched in " << q_nn_nj*1.0/process_unsuccessful*100 <<"%" << endl;
	cout<< "\t-neutralinos are on regular distance, their jets are regular, but more than 2 quarks lost their jets in "<< u << " and that is " 
		<< u*1.0/process_unsuccessful*100 <<"%" << endl;

}
int returnIndexSorted(double v_var,vector<double> *fatjet_var){	//returning index of current value in ascending sorted array
	sort(fatjet_var->begin(),fatjet_var->end() );
	for(int i=0;i<fatjet_var->size();i++)
		if(abs(v_var-  (*fatjet_var)[i] )<0.0000001)
			return i;
	return -1;
}
int returnIndex(double v_var,vector<double> *fatjet_var){	
	for(int i=0;i<fatjet_var->size();i++)
		if(abs(v_var -  (*fatjet_var)[i]) <0.0000001)
			return i;
	return -1;
} 

void checkingIfJetsHaveHighestPt(vector<TLorentzVector> min_fatjet,vector<double> *fatjet_pt,int& jets_pt_are_correlated){
	//if both fatjets are on highest 2 spots, we increase jets_pt_are_correlated variable.
	if( (returnIndexSorted(min_fatjet[0].Pt(),fatjet_pt) == fatjet_pt->size()-1 && returnIndexSorted(min_fatjet[1].Pt(),fatjet_pt) ==fatjet_pt->size()-2 )
		|| (returnIndexSorted(min_fatjet[0].Pt(),fatjet_pt) ==fatjet_pt->size()-2 && returnIndexSorted(min_fatjet[1].Pt(),fatjet_pt) 
		== fatjet_pt->size()-1 )
		) jets_pt_are_correlated++;

}

void checkingIfJetsHaveHighestM(vector<TLorentzVector> min_fatjet,vector<TLorentzVector> fatjet_rest,vector<double> *fatjet_pt,int& jets_M_are_correlated, vector<float>& matched_mass, vector<float> nmatched_mass){
	//if both fatjets are on highest 2 spots, we increase jets_M_are_correlated variable.
	vector<float> vector_tmp;
	matched_mass.push_back(min_fatjet[0].M());
	matched_mass.push_back(min_fatjet[1].M());
	for( int i=0 ; i< fatjet_rest.size(); i++){
		nmatched_mass.push_back(fatjet_rest[i].M());
		vector_tmp.push_back(fatjet_rest[i].M());
	}
	vector_tmp.push_back(min_fatjet[0].M());
	vector_tmp.push_back(min_fatjet[1].M());
	sort(vector_tmp.begin(),vector_tmp.end() );
	if( (matched_mass[0] == vector_tmp[vector_tmp.size()-1] && matched_mass[1] == vector_tmp[vector_tmp.size()-2]) ||  (matched_mass[0] == vector_tmp[vector_tmp.size()-2] && matched_mass[1] == vector_tmp[vector_tmp.size()-1]) ) jets_M_are_correlated++;

	// ovo je ideja za proveru Mase, da li su oba mecovana jeta zapravona ona dva sa najvecom masom

} 


void checkingIfJetsHaveHighestMax(int index0,int index1,vector<float> *fatjet_var,int& jets_pt_are_correlated){
	
	vector<double> vector_tmp;
	double v0 = (*fatjet_var)[index0],v1= (*fatjet_var)[index1];
        //copy(fatjet_var->begin(), fatjet_var->end(), back_inserter(*vector_tmp ));
	
	for(int i=0;i<fatjet_var->size();i++) {
		vector_tmp.push_back( ((*fatjet_var)[i])   );
		
	}
	//for(int i=0;i<fatjet_var->size();i++) {
	//	cout<< (vector_tmp[i]) << " ";
	//}
	//cout << endl;
	//for(int i=0;i<fatjet_var->size();i++) {//cout << ((*fatjet_var)[i]) << endl;
		//vector_tmp->push_back( ((*fatjet_var)[i])   );
		//(*vector_tmp)[i]= (*fatjet_var)[i];
	//	cout<< (vector_tmp[i]) << " ";
	//}
	//cout << endl;
	//cout << index0 << ", " << v0 << ", " << index1 << ", " << v1 << endl;
	sort(vector_tmp.begin(),vector_tmp.end() );	
	
	for(int i=0; i<vector_tmp.size();i++){
		if( abs( v0 - (vector_tmp[i]) )  <0.0000000001 ) index0 = i;		
		if( abs( v1 - (vector_tmp[i]) )  <0.0000000001 ) index1 = i;
	}
	//for(int i=0;i<fatjet_var->size();i++) {cout<< (vector_tmp[i]) << " ";
	//}
	//cout << index0 << ", " << v0 << ", " << index1 << ", " << v1 << endl;
	//cout << "-------------------------" << endl;
	/*
	if( (returnIndexSorted(v0,&vector_tmp) == vector_tmp.size()-1 && returnIndexSorted(v1,&vector_tmp) == vector_tmp.size()-2 )
		|| (returnIndexSorted(v0,&vector_tmp) == vector_tmp.size()-2 && returnIndexSorted(v1,&vector_tmp) == vector_tmp.size()-1)
		) jets_pt_are_correlated++;	*/
	if( ( index0 == vector_tmp.size()-1 && index1 == vector_tmp.size()-2 )
		|| ( index0 == vector_tmp.size()-2 && index1 == vector_tmp.size()-1)
		|| ( index0 == vector_tmp.size()-1 && index1 == vector_tmp.size()-1)
		 
		) jets_pt_are_correlated++;
 

} 
void checkingIfJetsHaveHighestMaxPt(int index0,int index1,vector<double> *fatjet_var,int& jets_pt_are_correlated){
	
	vector<double> vector_tmp;
	double v0 = (*fatjet_var)[index0],v1= (*fatjet_var)[index1];
        //copy(fatjet_var->begin(), fatjet_var->end(), back_inserter(*vector_tmp ));
	
	for(int i=0;i<fatjet_var->size();i++) {
		vector_tmp.push_back( ((*fatjet_var)[i])   );
		
	}
	//for(int i=0;i<fatjet_var->size();i++) {
	//	cout<< (vector_tmp[i]) << " ";
	//}
	//cout << endl;
	//for(int i=0;i<fatjet_var->size();i++) {//cout << ((*fatjet_var)[i]) << endl;
		//vector_tmp->push_back( ((*fatjet_var)[i])   );
		//(*vector_tmp)[i]= (*fatjet_var)[i];
	//	cout<< (vector_tmp[i]) << " ";
	//}
	//cout << endl;
	//cout << index0 << ", " << v0 << ", " << index1 << ", " << v1 << endl;
	sort(vector_tmp.begin(),vector_tmp.end() );	
	
	for(int i=0; i<vector_tmp.size();i++){
		if( abs( v0 - (vector_tmp[i]) )  <0.0000000001 ) index0 = i;		
		if( abs( v1 - (vector_tmp[i]) )  <0.0000000001 ) index1 = i;
	}
	//for(int i=0;i<fatjet_var->size();i++) {cout<< (vector_tmp[i]) << " ";
	//}
	//cout << index0 << ", " << v0 << ", " << index1 << ", " << v1 << endl;
	//cout << "-------------------------" << endl;
	/*
	if( (returnIndexSorted(v0,&vector_tmp) == vector_tmp.size()-1 && returnIndexSorted(v1,&vector_tmp) == vector_tmp.size()-2 )
		|| (returnIndexSorted(v0,&vector_tmp) == vector_tmp.size()-2 && returnIndexSorted(v1,&vector_tmp) == vector_tmp.size()-1)
		) jets_pt_are_correlated++;	*/
	if( ( index0 == vector_tmp.size()-1 && index1 == vector_tmp.size()-2 )
		|| ( index0 == vector_tmp.size()-2 && index1 == vector_tmp.size()-1)
		|| ( index0 == vector_tmp.size()-1 && index1 == vector_tmp.size()-1)
		 
		) jets_pt_are_correlated++;
 

} 
void checkingIfJetsHaveHighestMin(int index0,int index1,vector<float> *fatjet_var,int& jets_pt_are_correlated){
	
	vector<double> vector_tmp;
	double v0 = (*fatjet_var)[index0],v1= (*fatjet_var)[index1];
        //copy(fatjet_var->begin(), fatjet_var->end(), back_inserter(*vector_tmp ));
	
	for(int i=0;i<fatjet_var->size();i++) {
		vector_tmp.push_back( ((*fatjet_var)[i])   );
		
	}
	//for(int i=0;i<fatjet_var->size();i++) {
	//	cout<< (vector_tmp[i]) << " ";
	//}
	//cout << endl;
	//for(int i=0;i<fatjet_var->size();i++) {//cout << ((*fatjet_var)[i]) << endl;
		//vector_tmp->push_back( ((*fatjet_var)[i])   );
		//(*vector_tmp)[i]= (*fatjet_var)[i];
	//	cout<< (vector_tmp[i]) << " ";
	//}
	//cout << endl;
	//cout << index0 << ", " << v0 << ", " << index1 << ", " << v1 << endl;
	sort(vector_tmp.begin(),vector_tmp.end() );	
	
	for(int i=0; i<vector_tmp.size();i++){
		if( abs( v0 - (vector_tmp[i]) )  <0.0000000001 ) index0 = i;		
		if( abs( v1 - (vector_tmp[i]) )  <0.0000000001 ) index1 = i;
	}
	//for(int i=0;i<fatjet_var->size();i++) {cout<< (vector_tmp[i]) << " ";
	////}
	//cout << index0 << ", " << v0 << ", " << index1 << ", " << v1 << endl;
	//cout << "-------------------------" << endl;
	
	if( (returnIndexSorted(v0,&vector_tmp) == 0 && returnIndexSorted(v1,&vector_tmp) == 1 )
		|| (returnIndexSorted(v0,&vector_tmp) == 1 && returnIndexSorted(v1,&vector_tmp) == 0)
		|| (returnIndexSorted(v0,&vector_tmp) == 0 && returnIndexSorted(v1,&vector_tmp) == 0)

		) jets_pt_are_correlated++;	
	/*
	if( ( index0 == 0 && index1 == 1)
		|| ( index0 == 1 && index1 == 0)
		) jets_pt_are_correlated++; */

}

void jetCorrelationPrint(int jets_pt_are_correlated,int nentries,int process_unsuccessful){
	cout << "Number of proceses where jets matched with neutralino have highest pt:"  << jets_pt_are_correlated << endl;
	cout << "That makes " <<jets_pt_are_correlated*100.0/(nentries-process_unsuccessful) << "%"<< endl;
}

void overlap_removal(vector<TLorentzVector>& jet_pass, vector<TLorentzVector> min_fatjet,vector<double> *pt,vector<double> *eta,vector<double> *phi,vector<double> *e ){
	TLorentzVector tmp;
	for(int i=0; i< pt->size();i++){
		tmp.SetPtEtaPhiE((*pt)[i],(*eta)[i],(*phi)[i],(*e)[i]);
		if( tmp.DeltaR(min_fatjet[0]) > 1 || tmp.DeltaR(min_fatjet[1]) > 1 ){
			jet_pass.push_back(tmp);
		}
	}
}



void canvasDrawAndPrint(TH1F* histo_s,TH1F* histo_u,  TH1F* histo_tau1,TH1F* histo_tau2,TH1F* histo_tau3,TH1F* histo_tau21,TH1F* histo_tau32,TH1F* histo_tau31,TH1F* histo_tau1n,TH1F* histo_tau2n,TH1F* histo_tau3n,TH1F* histo_tau21n,TH1F* histo_tau32n,TH1F* histo_tau31n, TH1F* histo_m, TH1F* histo_mn){ //canvas draw
	TLegend *legendtau1,*legendtau2, *legendtau3, *legendtau21, *legendtau32, *legendtau31, *legendpt,*legendm;
	auto canvastau= new TCanvas();
	auto canvaspt= new TCanvas();
	auto canvasm= new TCanvas();

	canvastau->Divide(3,2);	
	canvastau->cd(1);
	histo_tau1n->SetLineColor(kRed);
   	histo_tau1n->DrawNormalized("Same");
	histo_tau1->SetLineColor(kBlue);
   	histo_tau1->DrawNormalized("Same");
	legendtau1= new TLegend(.7,.7,.9,.9,"successful events");
	
	legendtau1->SetFillColor(0);
	legendtau1->AddEntry(histo_tau1n,"Not matched fatjets");
	legendtau1->AddEntry(histo_tau1,"Matched fatjets");
	legendtau1->Draw("Same");
	
	canvastau->cd(2);
	histo_tau2n->SetLineColor(kRed);
   	histo_tau2n->DrawNormalized("Same");
   	histo_tau2->SetLineColor(kBlue);
	histo_tau2->DrawNormalized("Same");
	
	legendtau2= new TLegend(.7,.7,.9,.9,"successful events");
	
	legendtau2->SetFillColor(0);
	legendtau2->AddEntry(histo_tau2n,"Not matched fatjets");
	legendtau2->AddEntry(histo_tau2,"Matched fatjets");
	legendtau2->Draw("Same");

	canvastau->cd(3);
   	histo_tau3->SetLineColor(kBlue);
	histo_tau3->DrawNormalized("Same");
	histo_tau3n->SetLineColor(kRed);
   	histo_tau3n->DrawNormalized("Same");
	legendtau3= new TLegend(.7,.7,.9,.9,"successful events");

	legendtau3->SetFillColor(0);
	legendtau3->AddEntry(histo_tau3n,"Not matched fatjets");
	legendtau3->AddEntry(histo_tau3,"Matched fatjets");
	legendtau3->Draw("Same");

	canvastau->cd(4);
	histo_tau21n->SetLineColor(kRed);
   	histo_tau21n->DrawNormalized("Same");
   	histo_tau21->SetLineColor(kBlue);
	histo_tau21->DrawNormalized("Same");
	legendtau21= new TLegend(.7,.7,.9,.9,"successful events");
	
	legendtau21->SetFillColor(0);
	legendtau21->AddEntry(histo_tau21n,"Not matched fatjets");
	legendtau21->AddEntry(histo_tau21,"Matched fatjets");
	legendtau21->Draw("Same");

	canvastau->cd(5);
	histo_tau32n->SetLineColor(kRed);
   	histo_tau32n->DrawNormalized("Same");
   	histo_tau32->SetLineColor(kBlue);
	histo_tau32->DrawNormalized("Same");
	legendtau32= new TLegend(.7,.7,.9,.9,"successful events");
	
	legendtau32->SetFillColor(0);
	legendtau32->AddEntry(histo_tau32n,"Not matched fatjets");
	legendtau32->AddEntry(histo_tau32,"Matched fatjets");
	legendtau32->Draw("Same");

	canvastau->cd(6);
	histo_tau31->SetLineColor(kBlue);
	histo_tau31->DrawNormalized("Same");
	histo_tau31n->SetLineColor(kRed);
   	histo_tau31n->DrawNormalized("Same");
	
	
	legendtau31= new TLegend(.7,.7,.9,.9,"successful events");
	
	legendtau31->SetFillColor(0);
	legendtau31->AddEntry(histo_tau31n,"Not matched fatjets");
	legendtau31->AddEntry(histo_tau31,"Matched fatjets");
	legendtau31->Draw("Same");

	canvaspt->cd();
        histo_u->SetLineColor(kRed);
	histo_u-> DrawNormalized("Same");
        histo_s->SetLineColor(kBlue);
	histo_s-> DrawNormalized("Same");
	legendpt= new TLegend(.7,.7,.9,.9,"successful events");
	
	legendpt->SetFillColor(0);
	legendpt->AddEntry(histo_u,"Not matched fatjets");
	legendpt->AddEntry(histo_s,"Matched fatjets");
	legendpt->Draw("Same");

	canvasm->cd();
        histo_mn->SetLineColor(kRed);
	histo_mn-> DrawNormalized("Same");
        histo_m->SetLineColor(kBlue);
	histo_m-> DrawNormalized("Same");
	legendm= new TLegend(.7,.7,.9,.9,"successful events");
	
	legendm->SetFillColor(0);
	legendm->AddEntry(histo_mn,"Not matched fatjets");
	legendm->AddEntry(histo_m,"Matched fatjets");
	legendm->Draw("Same");

	//canvastau->Print("tau_histo1_1_nozero.pdf");
	//canvaspt->Print("pt_histo.pdf");
	//canvasm->Print("fatjet_mass.pdf");
	
}
void tauPrint(int process_successful, int jets_tau1_are_correlated=0, int jets_tau2_are_correlated=0, int jets_tau3_are_correlated=0, 
		int jets_tau21_are_correlated=0, int jets_tau31_are_correlated=0, int jets_tau32_are_correlated=0,
		int jets_NTrimSubjets_are_correlated=0,int jets_ungrtrk500_are_correlated=0){
	cout << "Tau1 correlation: " << (100.0*jets_tau1_are_correlated)/process_successful <<"%  #"<<jets_tau1_are_correlated<<endl;
	cout << "Tau2 correlation: " << (100.0*jets_tau2_are_correlated)/process_successful <<"%  #"<<jets_tau2_are_correlated<<endl;
	cout << "Tau3 correlation: " << (100.0*jets_tau3_are_correlated)/process_successful <<"%  #"<<jets_tau3_are_correlated<<endl;
	cout << "Tau21 correlation: " << (100.0*jets_tau21_are_correlated)/process_successful <<"%  #"<<jets_tau21_are_correlated<<endl;
	cout << "Tau31 correlation: " << (100.0*jets_tau31_are_correlated)/process_successful <<"%  #"<<jets_tau31_are_correlated<<endl;
	cout << "Tau32 correlation: " << (100.0*jets_tau32_are_correlated)/process_successful <<"%  #"<<jets_tau32_are_correlated<<endl;
	cout << "NTrimSubjets correlation: " << (100.0*jets_NTrimSubjets_are_correlated)/process_successful <<"%  #"<<jets_NTrimSubjets_are_correlated<<endl;
	cout << "ungrtrk500 correlation: " << (100.0*jets_ungrtrk500_are_correlated)/process_successful <<"%  #"<<jets_ungrtrk500_are_correlated<<endl;



}
void NTrimSubjets_ungrtrk500Print(TH1F *histo_NTrimSubjets,TH1F *histo_NTrimSubjetsn,TH1F *histo_ungrtrk500,TH1F *histo_ungrtrk500n){
	TLegend *legend_NTrimSubjets,*legend_ungrtrk500;
	auto canvas= new TCanvas();


	canvas->Divide(2,1);	
	canvas->cd(1);
	histo_NTrimSubjets->SetLineColor(kBlue);
   	histo_NTrimSubjets->DrawNormalized("Same");
	histo_NTrimSubjetsn->SetLineColor(kRed);
   	histo_NTrimSubjetsn->DrawNormalized("Same");
	legend_NTrimSubjets= new TLegend(.7,.7,.9,.9,"successful events");
	
	legend_NTrimSubjets->SetFillColor(0);
	legend_NTrimSubjets->AddEntry(histo_NTrimSubjetsn,"Not matched fatjets");
	legend_NTrimSubjets->AddEntry(histo_NTrimSubjets,"Matched fatjets");
	legend_NTrimSubjets->Draw("Same");
	
	canvas->cd(2);
	histo_ungrtrk500->SetLineColor(kBlue);
   	histo_ungrtrk500->DrawNormalized("Same");
	histo_ungrtrk500n->SetLineColor(kRed);
   	histo_ungrtrk500n->DrawNormalized("Same");
	legend_ungrtrk500= new TLegend(.7,.7,.9,.9,"successful events");
	
	legend_ungrtrk500->SetFillColor(0);
	legend_ungrtrk500->AddEntry(histo_ungrtrk500n,"Not matched fatjets");
	legend_ungrtrk500->AddEntry(histo_ungrtrk500,"Matched fatjets");
	legend_ungrtrk500->Draw("Same");

	//canvas->Print("NTRimSubjects_ungrtrk.pdf");
	

	
}
//MAIN
void Nova::Loop()
{	
//declaration   	
	if (fChain == 0) return;
	Long64_t nbytes = 0, nb = 0, nentries = fChain->GetEntriesFast();
	
   	vector<TLorentzVector> vector_neutralino;
	vector<TLorentzVector> min_fatjet(2);
	vector<TLorentzVector> fatjet_rest;
	vector<TLorentzVector> jet_pass;
	vector<TLorentzVector> quarks;	

   	TLorentzVector fatjet,quark;

   	int tmp_lost_quarks=0,total_lost_quarks=0;
   	int index=0,i=0,q=0,nn=0,nj=0,q_nn=0,q_nj=0,nn_nj=0,q_nn_nj=0,u=0;
   	int process_unsuccessful=0;
   	int jets_M_are_correlated=0,jets_pt_are_correlated=0,jets_tau1_are_correlated=0, jets_tau2_are_correlated=0, jets_tau3_are_correlated=0, jets_tau21_are_correlated=0, jets_tau31_are_correlated=0, jets_tau32_are_correlated=0;
	int jets_NTrimSubjets_are_correlated=0, jets_ungrtrk500_are_correlated=0;
  	bool process_unsuccessful_flag=false, provera=false;
	int tau1_zero_counter=0,tau2_zero_counter=0,tau3_zero_counter=0;
	float t=0;
	int good_match=0;
	int good_match_all=0;

	vector<float> *fatjet_tau31_wta, matched_mass, nmatched_mass;

	auto histo_u=new TH1F("not matched",
		"Fatjets matched/not matched with Neutralino; Pt[GeV] ; # of occurencies ",
		60, // Number of Bins
		0, // Lower X Boundary
		3000); // Upper X Boundary
	histo_u->SetStats(0);
	auto histo_s=new TH1F("mathced",
		"Fatjets matched/not matched with Neutralino; Pt[GeV] ; # of occurencies ",
		60, // Number of Bins
		0, // Lower X Boundary
		3000); // Upper X Boundary
	histo_s->SetStats(0);
	
	auto histo_tau1 = new TH1F("tau1",
			"Tau1;value;# of occurencies",
			30,
			0,
			1);
	histo_tau1->SetStats(0);
	auto histo_tau1n = new TH1F("tau1n",
			"Tau1;value;# of occurencies",
			30,
			0,
			1);
	histo_tau1n->SetStats(0);
	auto histo_tau2 = new TH1F("tau2",
			"Tau2;value;# of occurencies",
			30,
			0,
			1);
	histo_tau2->SetStats(0);
	auto histo_tau2n = new TH1F("tau2n",
			"Tau2;value;# of occurencies",
			30,
			0,
			1);
	histo_tau2n->SetStats(0);
	auto histo_tau3 = new TH1F("tau3",
			"Tau3;value;# of occurencies",
			30,
			0,
			1);
	histo_tau3->SetStats(0);
	auto histo_tau3n = new TH1F("tau3n",
			"Tau3;value;# of occurencies",
			30,
			0,
			1);
	histo_tau3n->SetStats(0);
	auto histo_tau21 = new TH1F("tau21",
			"Tau21;value;# of occurencies",
			30,
			0,
			1);
	histo_tau21->SetStats(0);
	auto histo_tau21n = new TH1F("tau21n",
			"Tau21;value;# of occurencies",
			30,
			0,
			1);
	histo_tau21n->SetStats(0);
	auto histo_tau32 = new TH1F("tau32",
			"Tau32;value;# of occurencies",
			30,
			0,
			1);
	histo_tau32->SetStats(0);
	auto histo_tau32n = new TH1F("tau32n",
			"Tau32;value;# of occurencies",
			30,
			0,
			1);
	histo_tau32n->SetStats(0);
	auto histo_tau31 = new TH1F("tau31",
			"Tau31;value;# of occurencies",
			30,
			0,
			1);
	histo_tau31->SetStats(0);
	auto histo_tau31n = new TH1F("tau31n",
			"Tau31;value;# of occurencies",
			30,
			0,
			1);
	histo_tau31n->SetStats(0);
    

	auto histo_NTrimSubjets = new TH1F("NTrimSubjets",
			"NTrimSubjets;value;# of occurencies",
			20,
			0,
			20);
	histo_NTrimSubjets->SetStats(0);
	auto histo_NTrimSubjetsn = new TH1F("NTrimSubjetsn",
			"NTrimSubjets;value;# of occurencies",
			20,
			0,
			20);
	histo_NTrimSubjetsn->SetStats(0);
 	auto histo_ungrtrk500 = new TH1F("ungrtrk500",
			"ungrtrk500;value;# of occurencies",
			150,
			0,
			150);
	histo_ungrtrk500->SetStats(0);
	auto histo_ungrtrk500n = new TH1F("ungrtrk500n",
			"ungrtrk500;value;# of occurencies",
			150,
			0,
			150);
	histo_ungrtrk500n->SetStats(0);
    
	auto histo_m = new TH1F("mass",
			"fatjet_mass;value[GeV];# of occurencies",
			120,
			-20,
			100);
	histo_m->SetStats(0);
	auto histo_mn = new TH1F("massn",
			"fatjet_mass;value[GeV];# of occurencies",
			120,
			-20,
			100);
	histo_mn->SetStats(0);
   
	SetActiveBranches(fChain);
	int grge=0;
//loop through processes
   	for (Long64_t jentry=0; jentry<nentries;jentry++){
        	Long64_t ientry = LoadTree(jentry);
        	if (ientry < 0) break;
        	nb = fChain->GetEntry(jentry);   nbytes += nb;
		//vector<int>*     truth_NeutralinoFromGluino_ParentBarcode
		provera=false;
		/*for( int i=0; i<fatjet_tau1_wta->size(); i++){
			if( (*fatjet_tau1_wta)[i]==t || (*fatjet_tau2_wta)[i]==t || (*fatjet_tau3_wta)[i]==t){
				provera=true;
			}
		}
		if (provera){
			continue;
		}*/

		overlap_removal(jet_pass,min_fatjet,jet_pt,jet_eta,jet_phi,jet_e );
		if(quark_from_gluino_matching(quarks,truth_QuarkFromGluino_pt,truth_QuarkFromGluino_eta,truth_QuarkFromGluino_phi,truth_QuarkFromGluino_e,jet_pass)==true){
			good_match++;
		}
		settingNeutralinos(vector_neutralino,truth_NeutralinoFromGluino_pt,
			truth_NeutralinoFromGluino_eta, truth_NeutralinoFromGluino_phi,truth_NeutralinoFromGluino_e);
		

        	findingMin(vector_neutralino,fatjet_pt, fatjet_eta,fatjet_phi,fatjet_e,min_fatjet,fatjet_rest);
		for(int i = 0;i<6;i++){
			(i<3)?(index=0):(index=1);
			quark.SetPtEtaPhiE((*truth_QuarkFromNeutralino_pt)[i],(*truth_QuarkFromNeutralino_eta)[i], 
				(*truth_QuarkFromNeutralino_phi)[i],(*truth_QuarkFromNeutralino_e)[i]);
		CheckIfQuarkIsOutOfJet(quark,min_fatjet[index],tmp_lost_quarks,total_lost_quarks,process_unsuccessful_flag);	
		}
		if(process_unsuccessful_flag==true){
			process_unsuccessful++;
			checkingForReasonOfProcessFailure(q,nn,nj,q_nn,q_nj,nn_nj,q_nn_nj,u,tmp_lost_quarks,vector_neutralino,min_fatjet);
			

		}else{
			if(quark_from_gluino_matching(quarks,truth_QuarkFromGluino_pt,truth_QuarkFromGluino_eta,truth_QuarkFromGluino_phi,truth_QuarkFromGluino_e,jet_pass)==true){
			good_match_all++;
			}
			int index0 =returnIndex(min_fatjet[0].Pt(),fatjet_pt);
			int index1 =returnIndex(min_fatjet[1].Pt(),fatjet_pt);
			//cout << min_fatjet[0].Pt() << " "<< (*fatjet_pt)[index0] << " " << (*fatjet_tau1_wta)[index0] << endl;
			
			histo_s->Fill(min_fatjet[0].Pt() );	//filling successful process histo
			histo_s->Fill(min_fatjet[1].Pt() );
			
			if((*fatjet_tau1_wta)[index0] !=t){
			histo_tau1->Fill( (*fatjet_tau1_wta)[index0] );
			}
			if((*fatjet_tau1_wta)[index1] !=t){  	//filling successful process histo
			histo_tau1->Fill( (*fatjet_tau1_wta)[index1] );
			}
			if((*fatjet_tau2_wta)[index0] !=t){
			histo_tau2->Fill((*fatjet_tau2_wta)[index0] );
			}	//filling successful process histo
			if((*fatjet_tau2_wta)[index1] !=t){
			histo_tau2->Fill((*fatjet_tau2_wta)[index1] );
			}
			if((*fatjet_tau3_wta)[index0] !=t){ 
			histo_tau3->Fill((*fatjet_tau3_wta)[index0] );	//filling successful process histo
			}
			if((*fatjet_tau3_wta)[index1] !=t){
			histo_tau3->Fill((*fatjet_tau3_wta)[index1] );
			}
			if((*fatjet_tau21_wta)[index0] !=t /*&& (*fatjet_tau21_wta)[index0] != "nan"*/ ){
			histo_tau21->Fill((*fatjet_tau21_wta)[index0]);	//filling successful process histo
			}
			if((*fatjet_tau21_wta)[index1] !=t /*&& (*fatjet_tau21_wta)[index1] != nan */){
			histo_tau21->Fill((*fatjet_tau21_wta)[index1]);	//filling successful process histo
			}
			if((*fatjet_tau32_wta)[index0] !=t /*&& (*fatjet_tau32_wta)[index0] != nan */ ){
			histo_tau32->Fill((*fatjet_tau32_wta)[index0]);	//filling successful process histo
			}
			if((*fatjet_tau32_wta)[index1] !=t /* && (*fatjet_tau32_wta)[index1] != nan */ ){
			histo_tau32->Fill((*fatjet_tau32_wta)[index1]);	//filling successful process histo
			}
					   //min_fatjet[0].Tau3()
			if(  ((*fatjet_tau3_wta)[index0]/(*fatjet_tau1_wta)[index0]) !=t /*&& ((*fatjet_tau3_wta)[index0]/(*fatjet_tau1_wta)[index0]) != nan */ ){
			histo_tau31->Fill((*fatjet_tau3_wta)[index0]/(*fatjet_tau1_wta)[index0] );	//filling successful process histo
			}
			if(  ((*fatjet_tau3_wta)[index1]/(*fatjet_tau1_wta)[index1]) !=t /* && ((*fatjet_tau3_wta)[index1]/(*fatjet_tau1_wta)[index1]) != nan */){
			histo_tau31->Fill((*fatjet_tau3_wta)[index1]/(*fatjet_tau1_wta)[index1] );
			}
			
			histo_NTrimSubjets->Fill((*fatjet_NTrimSubjets)[index0] );	//filling successful process histo
			histo_NTrimSubjets->Fill((*fatjet_NTrimSubjets)[index1] );
			
			histo_ungrtrk500->Fill((*fatjet_ungrtrk500)[index0] );	//filling successful process histo
			histo_ungrtrk500->Fill((*fatjet_ungrtrk500)[index1] ); 
			
			histo_m->Fill( min_fatjet[0].M() );
			histo_m->Fill( min_fatjet[1].M() );
			/*cout<< "min_fatjet[0].M()"<<" "<< min_fatjet[0].M()<<endl;
			cout<< "min_fatjet[1].M()"<<" "<< min_fatjet[1].M()<<endl;
			cout<< "0.Pt"<<" "<<"0.Eta"<<" "<< "0.Phi"<<" "<<"0.E"<<endl;
			cout<< (min_fatjet[0]).Pt()<<" "<<(min_fatjet[0]).Eta()<<" "<< (min_fatjet[0]).Phi()<<" "<< (min_fatjet[0]).E()<<endl;
			cout<< "1.Pt"<<" "<<"1.Eta"<<" "<< "1.Phi"<<" "<<"1.E"<<endl;
			cout<< min_fatjet[1].Pt()<<" "<<min_fatjet[1].Eta()<<" "<< min_fatjet[1].Phi()<<" "<<min_fatjet[1].E()<<endl;
			//cout<< "min_fatjet[1]"<<" "<< min_fatjet[1]<<endl;
			cout<<"--------------------------"<<endl;*/
			
			

			//for (int i =0; i<fatjet_pt->size();i++)
			//	cout << (*fatjet_pt)[i] << " ";
			//cout << endl;
//	fChain->SetBranchStatus("fatjet_ungrtrk500",1);
//	fChain->SetBranchStatus("fatjet_NTrimSubjets",1);

			checkingIfJetsHaveHighestMaxPt(index0,index1,fatjet_pt,jets_pt_are_correlated);
			
			checkingIfJetsHaveHighestMin(index0,index1,fatjet_tau1_wta,jets_tau1_are_correlated);
			checkingIfJetsHaveHighestMin(index0,index1,fatjet_tau2_wta,jets_tau2_are_correlated);
			checkingIfJetsHaveHighestMin(index0,index1,fatjet_tau3_wta,jets_tau3_are_correlated);
			checkingIfJetsHaveHighestMin(index0,index1,fatjet_tau21_wta,jets_tau21_are_correlated);
			fatjet_tau31_wta=new vector<float>();
			for(int i=0;i<fatjet_tau1_wta->size();i++)
				fatjet_tau31_wta->push_back( ((*fatjet_tau3_wta)[i]) / ((*fatjet_tau1_wta)[i])  );
			checkingIfJetsHaveHighestMin(index0,index1,fatjet_tau31_wta,jets_tau31_are_correlated);
			checkingIfJetsHaveHighestMin(index0,index1,fatjet_tau32_wta,jets_tau32_are_correlated);
			
			checkingIfJetsHaveHighestMax(index0,index1,fatjet_NTrimSubjets,jets_NTrimSubjets_are_correlated);
			checkingIfJetsHaveHighestMax(index0,index1,fatjet_ungrtrk500,jets_ungrtrk500_are_correlated);



/*			index0 =returnIndex(min_fatjet[0].Pt(),fatjet_pt);
			index1 =returnIndex(min_fatjet[1].Pt(),fatjet_pt);
			for (int i =0; i<fatjet_pt->size();i++)
				cout << (*fatjet_pt)[i] << " ";
			cout << endl;
			cout << "----------------------" << endl;
			sort(fatjet_pt->begin(),fatjet_pt->end() );
			for(int i=0;i<fatjet_pt->size();i++){
				if(abs(min_fatjet[0].Pt() - (*fatjet_pt)[i]) > 0.1){
					if(abs(min_fatjet[1].Pt() - (*fatjet_pt)[i]) > 0.1){
				//cout<< (*fatjet_pt)[i] <<" "<< min_fatjet[0].Pt()<<" " <<min_fatjet[1].Pt()<<endl;
						histo_u->Fill((*fatjet_pt)[i]);	
					}
			 	}
			}*/

			for(int i=0;i<fatjet_pt->size();i++){
				if(i!= index0 && i!= index1){
					histo_u->Fill((*fatjet_pt)[i]);	
					if( (*fatjet_tau1_wta)[i] !=t){
					histo_tau1n->Fill( (*fatjet_tau1_wta)[i] );	//filling unsuccessful process histo
					}
					if( (*fatjet_tau2_wta)[i] !=t){
					histo_tau2n->Fill( (*fatjet_tau2_wta)[i] );	//filling unsuccessful process histo
					}
					if( (*fatjet_tau3_wta)[i] !=t){
					histo_tau3n->Fill( (*fatjet_tau3_wta)[i] );	//filling unsuccessful process histo
					}
					if((*fatjet_tau21_wta)[i] !=t /*&& (*fatjet_tau21_wta)[i] != nan */){
					histo_tau21n->Fill((*fatjet_tau21_wta)[i]);	//filling unsuccessful process histo
					}
					if((*fatjet_tau32_wta)[i] !=t /*&& (*fatjet_tau32_wta)[i] != nan */){
					histo_tau32n->Fill((*fatjet_tau32_wta)[i]);	//filling unsuccessful process histo
					}
					if( (*fatjet_tau3_wta)[i]/(*fatjet_tau1_wta)[i]  !=t/* && (*fatjet_tau3_wta)[i]/(*fatjet_tau1_wta)[i]  != nan */){
					histo_tau31n->Fill((*fatjet_tau31_wta)[i]);	//filling unsuccessful process histo
					}
					//histo_tau31n->Fill((*fatjet_tau3_wta)[i]/(*fatjet_tau1_wta)[i] );
					histo_NTrimSubjetsn->Fill((*fatjet_NTrimSubjets)[i] );
	
					histo_ungrtrk500n->Fill((*fatjet_ungrtrk500)[i] );	//filling unsuccessful process histo
					if ( (*fatjet_ungrtrk500)[i] == 0 ) grge++;

				}
			}
			
			for(int i=0 ; i<fatjet_rest.size();i++){
				histo_mn->Fill( fatjet_rest[i].M() );
				//cout<< "fatjet_rest[i].M()"<<" "<< fatjet_rest[i].M()<<endl;
			//cout<< "min_fatjet[1].M()"<<" "<< min_fatjet[1].M()<<endl;
			//cout<< "Pt"<<" "<<"Eta"<<" "<< "Phi"<<" "<<"E"<<endl;
			//cout<< (fatjet_rest[i]).Pt()<<" "<<(fatjet_rest[i]).Eta()<<" "<< (fatjet_rest[i]).Phi()<<" "<< (fatjet_rest[i]).E()<<endl;
			//cout<< "1.Pt"<<" "<<"1.Eta"<<" "<< "1.Phi"<<" "<<"1.E"<<endl;
			//cout<< min_fatjet[1].Pt()<<" "<<min_fatjet[1].Eta()<<" "<< min_fatjet[1].Phi()<<" "<<min_fatjet[1].E()<<endl;
			//cout<< "min_fatjet[1]"<<" "<< min_fatjet[1]<<endl;
			//cout<<"--------------------------"<<endl;
			}
			checkingIfJetsHaveHighestM(min_fatjet,fatjet_rest,fatjet_pt,jets_M_are_correlated, matched_mass,nmatched_mass);
		}
		for( int i=0 ; i < fatjet_pt->size(); i++){
				if ( (*fatjet_tau1_wta)[i]==t){
					tau1_zero_counter++;
				}
			}	
			for( int i=0 ; i < fatjet_pt->size(); i++){
				if ( (*fatjet_tau2_wta)[i]==t ){
					tau2_zero_counter++;
				}
			}				
			for( int i=0 ; i < fatjet_pt->size(); i++){
				if ( (*fatjet_tau3_wta)[i]==t ){
					tau3_zero_counter++;
				}
			}
		resetFlags(process_unsuccessful_flag,tmp_lost_quarks,min_fatjet);
		fatjet_rest.clear();
		matched_mass.clear();
		nmatched_mass.clear();
		jet_pass.clear();
		quarks.clear();
		

		
	
   	}
//results 
	//cout<<  "tau1_zero_counter" << " "<< 100.0*tau1_zero_counter/(nentries)<<"%"<<endl;
	//cout<<  "tau2_zero_counter" << " "<< 100.0*tau2_zero_counter/(nentries)<<"%"<<endl;
	//cout<<  "tau3_zero_counter" << " "<< 100.0*tau3_zero_counter/(nentries)<<"%"<<endl;
	//cout<<"jets_M_are_correlated"<<" "<<jets_M_are_correlated<<" "<<"which is"<<" "<<100.0*jets_M_are_correlated/(nentries-process_unsuccessful)<<"%"<<endl;
	//NTrimSubjets_ungrtrk500Print(histo_NTrimSubjets,histo_NTrimSubjetsn,histo_ungrtrk500,histo_ungrtrk500n);
	//canvasDrawAndPrint(histo_s,histo_u, histo_tau1,histo_tau2,histo_tau3,histo_tau21,histo_tau32,histo_tau31, histo_tau1n,histo_tau2n,histo_tau3n,histo_tau21n,histo_tau32n,histo_tau31n,histo_m,histo_mn);	
	//jetCorrelationPrint(jets_pt_are_correlated,nentries,process_unsuccessful);
	//myPrint(total_lost_quarks,process_unsuccessful,nentries,nn,nj,q,nn_nj,q_nn,q_nj,q_nn_nj,u);
	//tauPrint(nentries-process_unsuccessful,jets_tau1_are_correlated,jets_tau2_are_correlated,jets_tau3_are_correlated,jets_tau21_are_correlated,jets_tau31_are_correlated,jets_tau32_are_correlated
	//	,jets_NTrimSubjets_are_correlated,jets_ungrtrk500_are_correlated);
	//cout <<"-----" << endl << grge << endl;
	cout<<"good_match"<<" "<<100.0*good_match/(nentries)<<endl;
	cout<<"good_match_all"<<" "<<100.0*good_match_all/(nentries)<<endl;
}







