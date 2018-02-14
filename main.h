void calculate_nonlinear_term(NonlinearTerm *bb, NonlinearTerm *vv, NonlinearTerm *bv, NonlinearTerm *vb,LinearTerm *b, LinearTerm *v, double DQx,double DQy,double DQz,double DVQ);

void calculate_linear_term(NonlinearTerm *bb, NonlinearTerm *vv, NonlinearTerm *bv, NonlinearTerm *vb, LinearTerm *b, LinearTerm *v, LinearTerm *fv, LinearTerm *fb, double DQx,double DQy,double DQz,double DVQ);

void integrate_over_t(LinearTerm *b, LinearTerm *v, LinearTerm *fv, LinearTerm *fb, double dt);

void set_initial_values_nonlinear_term(NonlinearTerm *nt);
void set_initial_values_linear_term(LinearTerm *lt,double DQx,double DQy,double DQz,double DVQ);
void set_initial_values_linear_term_gauss(LinearTerm *lt,double DQx,double DQy,double DQz,double DVQ);
void set_initial_values_nonlinear_term_zero(NonlinearTerm *nt);

