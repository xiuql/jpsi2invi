#include "PipiJpsiAlg/PionZeroList.h"
const double m_true_pi0 = 0.135;
using namespace std;

//====================
//********************
//====================
PionZero::PionZero(HepLorentzVector &vp_gam_1, int &index_gam_1, HepLorentzVector &vp_gam_2, int &index_gam_2):
  m_vp_1(vp_gam_1),
  m_index_1(index_gam_1),
  m_vp_2(vp_gam_2),
  m_index_2(index_gam_2){
  m_vp_total = m_vp_1 + m_vp_2;
  m_inv_mass = m_vp_total.m();
  cal_goodness(0);
}

////********************
//PionZero::PionZero():
//  m_vp_1(0),
//  m_index_1(x_gam_1),
//  m_vp_2(vp_gam_2),
//  m_index_2(index_gam_2){
//  m_vp_total = m_vp_1 + m_vp_2;
//  m_inv_mass = m_vp_total.m();
//  cal_goodness(0);
//}

//********************
void PionZero::cal_goodness( int good_method ){
  switch(good_method){
  case 0:
    m_goodness = fabs(m_inv_mass-m_true_pi0);
    break;
  case 1:
    cout << "not complete yet, sorry" << endl;
    break;
  default:
    cout << "this default do nothing" << endl;
  }
}

//====================
//********************
//====================
PionZeroList::PionZeroList(Vp4 &input_vp):  
  m_gam_vp(input_vp), 
  m_pi0_vp(0),
  Vpi0_list(0),
  m_low_cut(0.11),
  m_high_cut(0.15){
  int m_Ngam(input_vp.size());
  if(m_Ngam<2) cout << "number of photons less than 2" << endl;
  for(int i=0; i<m_Ngam-1; i++){
    for(int j=i+1; j<m_Ngam; j++){
      double m_temp_mass((m_gam_vp[i]+m_gam_vp[j]).m());
      if(m_temp_mass<m_low_cut || m_temp_mass>m_high_cut) continue;
      PionZero m_temp_pi0(m_gam_vp[i], i, m_gam_vp[j], j);
      Vpi0_list.push_back(m_temp_pi0);
    }
  }
  m_num_pi0 = Vpi0_list.size();
}

//*******************
void PionZeroList::set_cut(int index, double cut){
  if(index < 1) m_low_cut = cut;
  else m_high_cut = cut;
  std::vector<PionZero>::iterator m_ind;
  for(m_ind=Vpi0_list.begin(); m_ind!=Vpi0_list.end(); m_ind++){
    double m_temp_mass((*m_ind).get_mass());
    if(m_temp_mass<m_low_cut||m_temp_mass>m_high_cut) Vpi0_list.erase(m_ind); 
  }
  refresh();
}

//********************
void PionZeroList::sort(){ // sort the pion0 list respect to goodness
  std::sort(Vpi0_list.begin(), Vpi0_list.end(), f_less());
  refresh();
}

//********************
void PionZeroList::reduce(){
  std::vector<int> gam_index_v(0);
  std::vector<PionZero> temp_pi0_list(0);
  for(int i=0; i<Vpi0_list.size(); i++){
    int temp_one(Vpi0_list[i].get_index(0));
    int temp_two(Vpi0_list[i].get_index(1));
    if(find(gam_index_v.begin(), gam_index_v.end(), temp_one) != gam_index_v.end() ) continue;
    if(find(gam_index_v.begin(), gam_index_v.end(), temp_two) != gam_index_v.end() ) continue;
    gam_index_v.push_back(temp_one);
    gam_index_v.push_back(temp_two);
    temp_pi0_list.push_back(Vpi0_list[i]);
  }
  Vpi0_list = temp_pi0_list;
  refresh();
}

//********************
void PionZeroList::refresh(){
  m_num_pi0 = Vpi0_list.size();
  m_pi0_vp.clear();
  for(int i=0; i<m_num_pi0; i++) m_pi0_vp.push_back(Vpi0_list[i].get_pi0_vp());
}

//********************
void PionZeroList::print(){
  std::cout << "number of pion0 in the list " << m_num_pi0 << std::endl;
  for(int i=0; i<Vpi0_list.size(); i++){
    std::cout << "i= " << i << endl;
    std::cout << "pi0 mass " << Vpi0_list[i].get_mass() << endl; 
    std::cout << "pi0 goodness " << Vpi0_list[i].get_goodness() << endl;
    std::cout << "pi0 four momentum " << Vpi0_list[i].get_pi0_vp() << endl;
    std::cout << "fir gam index " << Vpi0_list[i].get_index(0) << endl;
    std::cout << "sec gam index " << Vpi0_list[i].get_index(1) << endl;
    std::cout << "fir gam mom " << Vpi0_list[i].get_gam_vp(0) << endl;
    std::cout << "sec gam mom " << Vpi0_list[i].get_gam_vp(1) << endl;
  }
}

