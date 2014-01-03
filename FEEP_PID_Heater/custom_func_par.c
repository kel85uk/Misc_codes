#include "udf.h"

#define NUM_UDM 3
static int udm_offset = UDM_UNRESERVED;
DEFINE_EXECUTE_ON_LOADING(on_loading, libname)
{
  if (udm_offset == UDM_UNRESERVED) udm_offset = Reserve_User_Memory_Vars(NUM_UDM);
  if (udm_offset == UDM_UNRESERVED)
    Message0("\nYou need to define up to %d extra UDMs in GUI and then reload current library %s\n", NUM_UDM, libname);
  else
  {
    Message0("%d UDMs have been reserved by the current library %s\n",NUM_UDM, libname);
    Set_User_Memory_Name(udm_offset,"Heater_power");
    Set_User_Memory_Name(udm_offset+1,"Error_PID");
    Set_User_Memory_Name(udm_offset+2,"lib1-UDM-2");
  }
  Message0("\nUDM Offset for Current Loaded Library = %d",udm_offset);
}

DEFINE_ON_DEMAND(set_udms)
{
  Domain *d;
  Thread *ct;
  cell_t c;
  int i;
  d=Get_Domain(1);
  if(udm_offset != UDM_UNRESERVED)
  {
    Message0("Setting UDMs\n");
    for (i=0;i<NUM_UDM;i++)
    {
      thread_loop_c(ct,d)
      {
        begin_c_loop(c,ct)
        {
          C_UDMI(c,ct,udm_offset+i)=0;
        }
        end_c_loop(c,ct)
      }
    }
  }
  else
    Message0("UDMs have not yet been reserved for library 1\n");
}

DEFINE_PROFILE(rad_elec_2_emit, t, i)
{
  Domain *domain;
  domain = Get_Domain(1);
  int ID=3; /* Get the ID for the extractor electrode */
  int axi = 1;
  Thread *tc=Lookup_Thread(domain,ID); /* Get the thread for cells in the extractor electrode */
  real vol_tot = 0;
  real temp_ave = 0;
  real axi_k = 0;
  cell_t c;  
  face_t f;
  
  axi_k = (axi==1)? 2*M_PI : 1;
  #if !RP_HOST
  begin_c_loop(c,tc) /* Obtain average temperature of the extractor electrode by volume averaging */
    {
      vol_tot += C_VOLUME(c,tc)*axi_k;
      temp_ave += C_T(c,tc)*C_VOLUME(c,tc)*axi_k;
    }
  end_c_loop(c,tc)
  #if RP_NODE
  temp_ave = PRF_GRSUM1(temp_ave);
  vol_tot = PRF_GRSUM1(vol_tot);
  temp_ave = temp_ave/vol_tot;
  
  Message0("Average extractor temperature = %f \n",temp_ave);

  begin_f_loop(f, t)
    {
      F_PROFILE(f, t, i) = temp_ave; /* assuming view factor from electrode to emitter surface = 1 */
    }
  end_f_loop(f, t)
  temp_ave = 0;
  vol_tot = 0;  
  #endif
  #if !PARALLEL
  temp_ave = temp_ave/vol_tot;
  Message0("Average extractor temperature = %f \n",temp_ave);

  begin_f_loop(f, t)
    {
      F_PROFILE(f, t, i) = temp_ave; /* assuming view factor from electrode to emitter surface = 1 */
    }
  end_f_loop(f, t)
  temp_ave = 0;
  vol_tot = 0;   
  #endif

  #endif
}

DEFINE_PROFILE(rad_elec_2_insul, t, i)
{
  Domain *domain;
  domain = Get_Domain(1);
  int ID=3; /* Get the ID for the extractor electrode */
  int axi = 1;
  Thread *tc=Lookup_Thread(domain,ID); /* Get the thread for cells in the extractor electrode */
  real vol_tot = 0;
  real temp_ave = 0;
  real axi_k = 0;
  cell_t c;  
  face_t f;
  
  axi_k = (axi==1)? 2*M_PI : 1;
  #if !RP_HOST
  begin_c_loop(c,tc) /* Obtain average temperature of the extractor electrode by volume averaging */
    {
      vol_tot += C_VOLUME(c,tc)*axi_k;
      temp_ave += C_T(c,tc)*C_VOLUME(c,tc)*axi_k;
    }
  end_c_loop(c,tc)
  #if RP_NODE
  temp_ave = PRF_GRSUM1(temp_ave);
  vol_tot = PRF_GRSUM1(vol_tot);
  temp_ave = temp_ave/vol_tot;
  Message0("Average extractor temperature = %f \n",temp_ave);

  begin_f_loop(f, t)
    {
      F_PROFILE(f, t, i) = temp_ave; /* assuming view factor from electrode to emitter surface = 1 */
    }
  end_f_loop(f, t)
  temp_ave = 0;
  vol_tot = 0;  
  #endif
  #if !PARALLEL
  temp_ave = temp_ave/vol_tot;
  Message0("Average extractor temperature = %f \n",temp_ave);

  begin_f_loop(f, t)
    {
      F_PROFILE(f, t, i) = temp_ave; /* assuming view factor from electrode to emitter surface = 1 */
    }
  end_f_loop(f, t)
  temp_ave = 0;
  vol_tot = 0;    
  #endif
  
  #endif
}

DEFINE_PROFILE(rad_emit_2_elec, t1, i1)
{
  Domain *domain1;
  domain1 = Get_Domain(1);
  int ID1=4; /* Get the ID for the insulator */
  int axi1 = 1;
  Thread *tc1=Lookup_Thread(domain1,ID1); /* Get the thread for cells in the insulator */
  real vol_tot1 = 0;
  real temp_ave1 = 0;
  real axi_k1 = 0;
  cell_t c1;  
  face_t f1;
  
  axi_k1 = (axi1==1)? 2*M_PI : 1;
  #if !RP_HOST
  begin_c_loop(c1,tc1) /* Obtain average temperature of the insulator by volume averaging */
    {
      vol_tot1 += C_VOLUME(c1,tc1)*axi_k1;
      temp_ave1 += C_T(c1,tc1)*C_VOLUME(c1,tc1)*axi_k1;
    }
  end_c_loop(c1,tc1)
  #if RP_NODE
  temp_ave1 = PRF_GRSUM1(temp_ave1);
  vol_tot1 = PRF_GRSUM1(vol_tot1);
  temp_ave1 = temp_ave1/vol_tot1;
  
  Message0("Average insulator temperature = %f \n",temp_ave1);

  begin_f_loop(f1, t1)
    {
      F_PROFILE(f1, t1, i1) = temp_ave1; /* assuming view factor from insulator to extractor surface = 1 */
    }
  end_f_loop(f1, t1)
  temp_ave1 = 0;
  vol_tot1 = 0;  
  #endif
  #if !PARALLEL
  temp_ave1 = temp_ave1/vol_tot1;
  Message0("Average insulator temperature = %f \n",temp_ave1);

  begin_f_loop(f1, t1)
    {
      F_PROFILE(f1, t1, i1) = temp_ave1; /* assuming view factor from insulator to extractor surface = 1 */
    }
  end_f_loop(f1, t1)
  temp_ave1 = 0;
  vol_tot1 = 0;   
  #endif
  
  #endif
}

DEFINE_SOURCE(heater_pow,c,t,dS,eqn) /* Constant heater power */
{
  Domain *domain1;
  domain1 = Get_Domain(1);
  int ID1=5; /* Get the ID for the heater */
  int axi1 = 1;
  Thread *tc1=Lookup_Thread(domain1,ID1); /* Get the thread for cells in the emitter */
  real vol_tot1 = 0;
  real heater_power = 0.5; /* in Watt */
  real axi_k1 = 0;
  real source = 0;
  cell_t c1;  
  face_t f1;
  
  axi_k1 = (axi1==1)? 2*M_PI : 1;
  
  #if !RP_HOST
  begin_c_loop(c1,tc1) /* Obtain average temperature of the extractor electrode by volume averaging */
    {
      vol_tot1 += C_VOLUME(c1,tc1)*axi_k1;
    }
  end_c_loop(c1,tc1)
  
  #if RP_NODE
  vol_tot1 = PRF_GRSUM1(vol_tot1);
  
/*  Message0("Volume of heater = %lf \n",vol_tot1); */
  source = heater_power/vol_tot1;
  dS[eqn] = 0;

  /* Enter into UDMI for post_processing */
  C_UDMI(c,t,0) = heater_power;
  return source;  
  #endif
  #if !PARALLEL
  vol_tot1 = vol_tot1;
/*  Message0("Volume of heater = %lf \n",vol_tot1); */
  source = heater_power/vol_tot1;
  dS[eqn] = 0;

  /* Enter into UDMI for post_processing */
  C_UDMI(c,t,0) = heater_power;
  return source;  
  #endif  

  #endif
}

DEFINE_SOURCE(heater_pow_P_SS,c,t,dS,eqn) /* Proportional-SS controller */
{
  Domain *domain1;
  domain1 = Get_Domain(1);
  int ID_heat=5; /* Get the ID for the heater */
  int ID_emit=6; /* Get the ID for the emitter */
  int axi1 = 1; /* 1 if axisymmetric, 0 if cartesian */
  Thread *tc_heat=Lookup_Thread(domain1,ID_heat); /* Get the thread for cells in the heater */
  Thread *tc_emit=Lookup_Thread(domain1,ID_emit); /* Get the thread for cells in the emitter */
  real vol_tot_heat = 0;
  real vol_tot_emit = 0;
  real temp_emit = 0;
  real temp_target = 473.15; /* 200 C */
  real heater_power_max = 0.5; /* in Watt */
  real heater_power_ss = 0;
  real dheat_dT = 0;
  real heater_power = 0;
  /* Controller constants */
  real kp = 0.1;
/*  real kd = 0.1;
  real ki = 0.1; */
  real axi_k1 = 0;
  real source = 0;
  cell_t c1;
  face_t f1;
  
  axi_k1 = (axi1==1)? 2*M_PI : 1;
  
  #if !RP_HOST
  heater_power_ss = 4.382e-6*pow(temp_target,2)-0.001093*temp_target-0.06257;

  begin_c_loop(c1,tc_heat) /* Obtain total heater volume */
    {
      vol_tot_heat += C_VOLUME(c1,tc_heat)*axi_k1;
    }
  end_c_loop(c1,tc_heat)
  
  begin_c_loop(c1,tc_emit) /* Obtain emitter temperature */
    {
      vol_tot_emit += C_VOLUME(c1,tc_emit)*axi_k1;
      temp_emit += C_T(c1,tc_emit)*C_VOLUME(c1,tc_emit)*axi_k1;
    }
  end_c_loop(c1,tc_emit)
  #if RP_NODE
  vol_tot_heat = PRF_GRSUM1(vol_tot_heat);
  vol_tot_emit = PRF_GRSUM1(vol_tot_emit);
  temp_emit = PRF_GRSUM1(temp_emit);
  temp_emit = temp_emit/vol_tot_emit;
  #endif
  #if !PARALLEL
  temp_emit = temp_emit/vol_tot_emit;
  #endif

  dheat_dT = 2*4.382e06*temp_emit-0.001093;
  kp = 2*4.382e06*temp_target-0.001093;
  heater_power = -kp*(temp_emit - temp_target) + heater_power_ss;
  heater_power = (heater_power >= 0)? heater_power : 0;
  heater_power = (heater_power<=heater_power_max)? heater_power: heater_power_max;
  source = heater_power/vol_tot_heat;
  dS[eqn] = 0;

  /* Enter into UDMI for post_processing */
  C_UDMI(c,t,0) = heater_power;
  #endif
  return source;
}

DEFINE_SOURCE(heater_pow_bang,c,t,dS,eqn) /* Bang-Bang controller minimum time */
{
  Domain *domain1;
  domain1 = Get_Domain(1);
  int ID_heat=5; /* Get the ID for the heater */
  int ID_emit=6; /* Get the ID for the emitter */
  int axi1 = 1; /* 1 if axisymmetric, 0 if cartesian */
  Thread *tc_heat=Lookup_Thread(domain1,ID_heat); /* Get the thread for cells in the heater */
  Thread *tc_emit=Lookup_Thread(domain1,ID_emit); /* Get the thread for cells in the emitter */
  real vol_tot_heat = 0;
  real vol_tot_emit = 0;
  real temp_emit = 0;
  real temp_target = 473.15; /* 200 C */
  real heater_power_max = 0.5;
  real heater_power = 0.5;
  real axi_k1 = 0;
  real source = 0;
  cell_t c1;
  face_t f1;
  
  axi_k1 = (axi1==1)? 2*M_PI : 1;
  #if !RP_HOST
  begin_c_loop(c1,tc_heat) /* Obtain total heater volume */
    {
      vol_tot_heat += C_VOLUME(c1,tc_heat)*axi_k1;
    }
  end_c_loop(c1,tc_heat)
  
  begin_c_loop(c1,tc_emit) /* Obtain emitter temperature */
    {
      vol_tot_emit += C_VOLUME(c1,tc_emit)*axi_k1;
      temp_emit += C_T(c1,tc_emit)*C_VOLUME(c1,tc_emit)*axi_k1;
    }
  end_c_loop(c1,tc_emit)
  #if RP_NODE
  vol_tot_heat = PRF_GRSUM1(vol_tot_heat);
  vol_tot_emit = PRF_GRSUM1(vol_tot_emit);
  temp_emit = PRF_GRSUM1(temp_emit);
  temp_emit = temp_emit/vol_tot_emit;
  #endif
  #if !PARALLEL
  temp_emit = temp_emit/vol_tot_emit;
  #endif
  heater_power = (temp_emit <= temp_target)? heater_power_max : 0;
  source = heater_power/vol_tot_heat;
  dS[eqn] = 0;
  /* Enter into UDMI for post_processing */
  C_UDMI(c,t,0) = heater_power;
  
  #endif
  return source;
}

DEFINE_SOURCE(heater_pow_PID,c,t,dS,eqn) /* PID controller */
{
  Domain *domain1;
  domain1 = Get_Domain(1);
  int ID_heat=5; /* Get the ID for the heater */
  int ID_emit=6; /* Get the ID for the emitter */
  int axi1 = 1; /* 1 if axisymmetric, 0 if cartesian */
  Thread *tc_heat=Lookup_Thread(domain1,ID_heat); /* Get the thread for cells in the heater */
  Thread *tc_emit=Lookup_Thread(domain1,ID_emit); /* Get the thread for cells in the emitter */
  real vol_tot_heat = 0;
  real vol_tot_emit = 0;
  real cell_vol = 0;
  real temp_emit = 0;
  real temp_emit_M1 = 0;
  real dtemp_dt = 0;
  real dt = CURRENT_TIMESTEP;
  real temp_target = 473.15; /* 200 C */
  real heater_power_max = 0.5; /* in Watt */
  real heater_power_ss = 0;
  real dheat_dT = 0;
  real heater_power = 0;
  /* Controller constants */
  real kp = 80*2.857e-3;
  real kd = 2*kp;
  real ki = 4e-5*kp;
  real axi_k1 = 0;
  real source = 0;
  cell_t c1;
  face_t f1;
  
  axi_k1 = (axi1==1)? 2*M_PI : 1;
  #if !RP_HOST
  begin_c_loop(c1,tc_heat) /* Obtain total heater volume */
    {
      vol_tot_heat += C_VOLUME(c1,tc_heat)*axi_k1;
    }
  end_c_loop(c1,tc_heat)
  
  begin_c_loop(c1,tc_emit) /* Obtain emitter temperature at current and previous timesteps */
    {
      cell_vol = C_VOLUME(c1,tc_emit)*axi_k1;
      vol_tot_emit += cell_vol;
      temp_emit += C_T(c1,tc_emit)*cell_vol;
      temp_emit_M1 += C_T_M1(c1,tc_emit)*cell_vol;
    }
  end_c_loop(c1,tc_emit)
    #if RP_NODE
  vol_tot_heat = PRF_GRSUM1(vol_tot_heat);
  vol_tot_emit = PRF_GRSUM1(vol_tot_emit);
  temp_emit = PRF_GRSUM1(temp_emit);
  temp_emit_M1 = PRF_GRSUM1(temp_emit);
  temp_emit = temp_emit/vol_tot_emit;
  temp_emit_M1 = temp_emit_M1/vol_tot_emit;
  dtemp_dt = (temp_emit - temp_emit_M1)/dt;
  #endif
  #if !PARALLEL
  temp_emit = temp_emit/vol_tot_emit;
  temp_emit_M1 = temp_emit_M1/vol_tot_emit;
  dtemp_dt = (temp_emit - temp_emit_M1)/dt;
  #endif

  heater_power = -kp*(temp_emit - temp_target) - kd*dtemp_dt + ki*C_UDMI(c,t,1);
  heater_power = (heater_power >= 0)? heater_power : 0;
  heater_power = (heater_power<=heater_power_max)? heater_power: heater_power_max;
  source = heater_power/vol_tot_heat;
  dS[eqn] = 0;
  /* Enter UDMI 1 for error integral of PID controller  C_UDMI(1) = int(e dt,0,t-1) */
  C_UDMI(c,t,1) = C_UDMI(c,t,1) + (temp_target-temp_emit)*dt;
  
  /* Enter into UDMI for post_processing */
  C_UDMI(c,t,0) = heater_power;
  #endif 
  return source;
}
