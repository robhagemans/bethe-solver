
#include "roots.h"
#include "nonsymmetric-string.h"
#include "solver.h"


#include "roots.h"
#include "symmetric-string.h"
#include "solver.h"

using namespace std;

void showIteration (StateSolver* solver)
{
    cout<< "iter "<< solver->iterations() <<" ";
    cout<<" conv "<< solver->convergence() <<" ";
    cout<<" rts "<< getValues(solver->getRoots()) <<" ";
    cout<<endl;
};



int main()
{
    double desired_convergence = 1e-20;
    
    int N;
    StateSolver* all_roots = 0;
    
    cout<<"test: ";
    int test;
    cin>>test;
    
    cout<<" maxiter: ";
    int max_iter = 200;
    cin>>max_iter;
    
    //for(int test=1; test<=100;++test) 
    {
        cout<<endl<<"test case #"<<test<<endl;
           
        switch (test) {
        case 1:
            //#0#1: 2I [ +1 -1 +3 +5 ;] bt roots (0.0485441,0) (-0.108553,0) (0.21904,0) (0.467018,0)
            //1;-1;3;5; 2J [+1 -1 +3 +5] sums [+1;-1;+3;+5;] sum +8 (+0.048544097,+0) (-0.10855296,+0) (+0.21903978,+0) (+0.46701764,+0)
           // qn (0.505131,0) (-0.491802,0) (1.50107,0) (2.48561,0)

            N = 12;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymRealRoots (N, {1, -1, 3, 5}) }, &showIteration);
            break;
        case 2:    
// state [1 2]#42
// 2I [[-6] [+3 -3]] 
//-10
//-6
//10
// iter 97
// conv 9.41352e-21
// roots [(-0.729231,0) (1.01393,0.52964) (1.01393,-0.52964) (-0.708918,0.506252) (-0.708918,-0.506252)]
// jx2 [-6 6 8 -8 -6]
// err 1.93067e-10

 
            // 2J [-6 +6 +8 -8 -6] sums [-6;-10;+10;] sum -6 (-0.72923137,+0) (+1.0139262,+0.52964034) (+1.0139262,-0.52964034) (-0.70891822,+0.50625222) (-0.70891822,-0.50625222)
            // roots (-0.729231369546914,0) (1.01392618798011,0.529640341317901) (1.01392618798011,-0.529640341317901) (-0.708918223625972,0.506252220106076) (-0.708918223625972,-0.506252220106076)

            //qn (-2.9999999999941,0) (3.00000000000053,-1.68984107108533e-11) (4.00000000000053,1.68984107108533e-11) (-4.00000000001264,-1.32917359150162e-10) (-3.00000000001264,1.32917359150162e-10)

            N=12;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymRealRoots (N, {-6}), new NonsymString(N, 2,-3), new NonsymString(N, 2,3) }, &showIteration);

//            all_roots = new Bunch<NonsymRoots> ({ new NonsymString (N, 1, -6), new NonsymString(N, 2,-3), new NonsymString(N, 2,3) }, &showIteration);

            break;
        case 3:
            N=10;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymRealRoots (N, {0, 2}), new NonsymString(N, 3, 0, 22) }, &showIteration);
            break;
        case 4:
            N=10;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymRealRoots (N, {0, 2}), new NonsymString(N, 3, 0, 2) }, &showIteration);
            break;
        case 5:
            N=10;
            //2I = [1, -5; 0 ]  2J = [2, -4; -2]
            //2J [+2 -4 +6 +4 +8] sums [+2;-4;-2;] sum -4 (+0.10226663,+0) (-0.74048306,+0) (+0.21276772,+1.0098234) (+0.21268098,+0) (+0.21276772,-1.0098234)

            //all_roots = { new NonsymString(1, 1, 2, 0.0791922) ,new NonsymString(1, -5, -4, -0.5), new NonsymString(3, 0, -2, 0.) };
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, 1) ,new NonsymString(N, 1, -5), new NonsymString(N, 3, 0) }, &showIteration);
            
            //all_roots = { new NonsymRealRoots ({2, -4}, {0.0791922, -0.5}), new NonsymString(3, 0, -2, 0.) };
            //all_roots = { new NonsymRealRoots ({2, -4}, {0.108377, -0.737097}), new NonsymString(3,-2, 0.200666) };
            
            break;
         case 6:
               //            #0#3: 2I [ +1 -3 ;; +0 ;] bt roots (0.101139,0) (-0.320041,0) (0.0835272,1) (0.0835272,0) (0.0835272,-1)
            // 2J [+0 -2 +8 +2 +10] sums [+0;-2;+0;] sum -2 (+0.12418901,+0) (-0.32285804,+0) (+0.072502535,+1.01871) (+0.053663952,+0) (+0.072502535,-1.01871)

            
            N=10; 
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, 1) ,new NonsymString(N, 1, -3), new NonsymString(N, 3, 0) }, &showIteration);
            /*
            NonsymString* str = new NonsymString(3, 0, 0,0.0835272);
            str->damping_delta_=0.0;
            all_roots = { new NonsymString(1, 1, 0, 0.101139) ,new NonsymString(1, -3, -2, -0.320041), str };
            */
            break;
         case 7:   
            // base=(1 1 1)  2I=(+4 -2 +0)  2J=(+3 +12 +23 )
// are those two roots exactly equal?
            
            N=12;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, 4) ,new NonsymString(N, 2, -2), new NonsymString(N, 3, 0) }, &showIteration);
            break;
         case 8:
            //   base=(1 1 0 1)  2I=(+2 +0  +0) 
            //guess-2J=(+0 +10 +40 )
            N=14;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, 2) ,new NonsymString(N, 2, 0), new NonsymString(N, 4, 0) }, &showIteration);
            break;
         case 9:
            N=14;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, 1) , new NonsymString(N, 1, -9), new NonsymString(N, 3, -4) }, &showIteration);
               
            //             base=(2 0 1)  2I=(+1 -9  -4) 
            break;
         case 10:
            N=12;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, -8), new NonsymString(N, 3, +4) }, &showIteration);
            break;
            // base=(1 0 1)  2I=(-8  +4)  2J=(-7 +27 )
            //conv 4.5687e-21
            //iter 30
            //roots (-1.3843,0) (1.19464,1.03026) (1.25182,0) (1.19464,-1.03026)
            //dirty (-3.5,0) (4.5,-3.7527e-11) (3.5,0) (5.5,3.75269e-11)
         case 11:
            N=12;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, -8), new NonsymString(N, 3, -4) }, &showIteration);
            break;
            //base=(1 0 1)  2I=(-8  -4)  2J=(-9 +21 )
            //conv 4.9385e-21
            //iter 46
            //roots (-0.756513,0) (-1.16252,0.978651) (-1.40483,0) (-1.16252,-0.978651)
            //dirty (-4.5,0) (-5.5,-2.98476e-11) (-4,0) (-4.5,2.98476e-11)
            // not a solution!
        
        
        // N=12, m=6 (0 -2 6  0) 
         case 12:
            N=12;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, 0), new NonsymString(N, 1, -2), new NonsymString(N, 1, 6), new NonsymString(N, 3, 0) }, &showIteration);
            break;
            // roots (-0.00481515,0) (-0.154716,0) (0.799635,0) (-0.207524,1.01958) (-0.225055,0) (-0.207524,-1.01958)

            //qn (-0.5,0) (-1.5,0) (2.5,0) (-5.5,-1.41402e-10) (-0.5,0) (-4.5,1.41402e-10)
            // converges with limited crossing (epsilons cannot move across deltas), no recalculation of Js.
            // also converges with no crossing        
            break;
         case 13:
            N=12; //M=4;
            //(3 7 4)
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, 3), new NonsymString(N, 1, 7), new NonsymString(N, 2, 4) }, &showIteration);
            // does not converge with no crossing: stalls
            //3 7 -8
            //iter 999  conv 2  rts (0.124642,0) (0.713826,0) (0.713826,0.509028) (0.713826,-0.509028)
            
            // does converge with 'limited crossing':
            // 3 7 -8
            // roots (0.124821,0) (0.719458,0) (0.708548,0.508552) (0.708548,-0.508552)

            //qn (1.5,0) (3.5,0) (3.5,9.90658e-11) (4.5,-9.90658e-11)
            // bethe-takahashi soln: iter 35  conv 2  rts (0.123231,0) (0.695901,0) (0.749503,0.5) (0.749503,-0.5) 
            // with j = 3 7 -8. how can j be the same if the order of strings is different?
            // because sign term is zero between 1 and 2-strings. 
            // => so 1-strings can pass 2-strings!
            break;
        case 14:
            N=12;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N,1, -6), new NonsymString(N, 3, -2) }, &showIteration);
            break;
        
        case 15:   
            /*
             state [1 1 1]#19
 2I [[-4] [-2] [+0]] 
-3
12
-3
 iter 49
 conv 6.56439e-21
 roots [(-0.303681,0) (-0.641278,0.50335) (-0.641278,-0.50335) (0.529579,0.996605) (0.527081,0) (0.529579,-0.996605)]
 jx2 [-3 -7 -5 7 5 9]
 err 1.57954e-10

            */
            
            N=12;
            all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, -4) ,new NonsymString(N, 2, -2), new NonsymString(N, 3, 0) }, &showIteration);
            
            
            break;

        }


        if (!all_roots) {
            cout<<" no test case with that number"<<endl;
            return 1;
        }
        
        IterResult solved = all_roots->solve(max_iter, desired_convergence);
        vector<CVal> roots = all_roots->getRoots();
        delete all_roots;
        
        
        vector<CVal> roots1 = roots;//all_roots->getRoots();
            vector< vector<CVal> > roots_string = {
                vector<CVal> { roots1[0] },
                vector<CVal> { roots1[1], roots1[2] },
                vector<CVal> { roots1[3], roots1[4], roots1[5] }
                };
            cout<<"I_BT? "<<getIs_BT(N, roots_string)<<endl;            
            
        
        cout.precision(15);
        cout<<endl<<" test "<<test<<endl;
        cout<<" N "<<N<<endl; 
        cout<<" roots "<<getValues(roots)<<endl;
        cout.precision(4);
        cout<<endl<<" dirty "<< getDirtyJs(N, roots) <<endl;
        cout<<endl<<" clean "<< getJs(N, roots) <<endl;
        cout<<endl<<" BT "<< getJs_BT(N, roots) <<endl;

        double betherr = 0;
        bool parityerr = false;
        cout<<endl<<" jx2 "<< getJx2(N, roots, &betherr, &parityerr) <<endl;
        cout<<" err "<< betherr<<endl;
        cout<<endl;
        if (solved != IterResult::iter_ok) 
            cout<<"\033[1;31m "<<message(solved)<<" \033[0m"<<endl; 
        if (parityerr) 
            cout<<"\033[1;31m parity error \033[0m"<<endl; 
        if (betherr>1e-5)
            cout<<"\033[1;31m large quantum number error \033[0m"<<endl;
    }
    return 0;    
}




