deviance(j,alpha,a) = -0.5;
				deviance(j,alpha,a+1) = 0.5;

				// bethe qns should be equal for collapsing string
				// why not I/2 ? ah well, I'm getting old I guess
				double bethe_quantum_x2 = -19;//quantum_number(j,alpha)/2;
cerr<<	bethe_quantum_x2 <<endl;
				double rap = rapidity(j,alpha);
				double eps = 0;
	cerr<<'*' <<SEP<< roots()<<endl;
				// hold the roots.
				hold(j, alpha) = true;

				for (int i=0; i<1000; ++i ) {

					// try for two-string
					double new_rap = 0.5* tan( 1.0/double(p_chain->length()) * ( 0.5*PI* bethe_quantum_x2 +  atan(-eps) ) );
					double new_eps = eps - p_chain->length() * atan(2.0*rap) - 0.5*PI* bethe_quantum_x2- atan(eps) ;
//cerr<<new_rap<<SEP<<

					for (int stepx = -100; stepx<100;++stepx) {
						double x = stepx/10.0;

						double last = -100;
						for (int stepy = -10; stepy<10;++stepy) {
							double y = rap + stepy/10.0;

							//cout<<x<< SEP<<y<<SEP<<
							if (sgn(last) != sgn(double(p_chain->length()) * atan(2.0*y) - ( 0.5*PI* bethe_quantum_x2 +  atan(-x) ))) {
								cout<<x<<SEP<<y<<last<<SEP<<endl;

							}
							last = double(p_chain->length()) * atan(2.0*y) - ( 0.5*PI* bethe_quantum_x2 +  atan(-x) );
						}
						cout<<endl;
					}
cerr<<bethe_quantum_x2<<endl;
					exit(0);

					rap = rapidity(j,alpha) = 0.1*rap + 0.9*new_rap;
					eps = aberration(j,alpha,a) = new_eps;
					//aberration(j,alpha,a) = 0;

					//iterate();
	cerr<<'*' <<SEP<< roots()<<endl;
				}