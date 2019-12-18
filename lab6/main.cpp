#include <iostream>
#include <vector>
#include <random>
#include <iomanip>

#include "function.h"

struct Answer
{
    double lambda;//==h
    double dis;//==final criterij
    std::vector <double> alpha;
    double w;
    double d;
    double J;


};
void disp1(double h, double J, double w, double d){
    std::cout << std::endl << std::endl << "    ------------------------------------------ " << std::endl;
    std::cout << std::fixed <<std::setprecision(4) <<std::setw(7)<< " || "<<
                            std::setw(6) << 'h' <<" | " <<
                            std::setw(7) << 'J' << " | " <<
                            std::setw(7) << 'w' << " | " <<
                            std::setw(7) << 'd' << " || "<< std::endl;
    std::cout << "    ------------------------------------------ " << std::endl;
    std::cout << std::fixed <<std::setprecision(4) <<std::setw(7)<< " || " <<
                            std::setw(2) << h << " | " <<
                            std::setw(7) << J << " | " <<
                            std::setw(7) << w << " | " <<
                            std::setw(7) << d << " ||"<< std::endl;
    std::cout << "    ------------------------------------------ " << std::endl << std::endl << std::endl;
}
void head(int a){
    std::cout << "    ";
    for(size_t i = 0; i < a; ++i) {
        std::cout << "----------";
    }
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << std::fixed <<std::setprecision(4) <<std::setw(7)<< " || "<<
              std::setw(6) << 'h' <<" | " <<
              std::setw(7) << "dis" << " | ";

    for(size_t i = 0; i < a; ++i){
        if(i==a/2) {
            std::cout << std::fixed << std::setprecision(4) <<
                      std::setw(7) << "   alpha  ";
        }
        else std::cout << std::fixed << std::setprecision(4) <<
                       std::setw(7) << "            ";
    }

    std::cout << std::fixed <<std::setprecision(4) <<
                        std::setw(7) << " |     w  " << " | " <<
                        std::setw(7) << "     d" << " || "<< std::endl<<"    ";
    for(size_t i = 0; i < a; ++i) {
        std::cout << "----------";
    }
    std::cout << "-------------------------------------------------" << std::endl;
}

void head2(int a){
    std::cout << "    ";
    for(size_t i = 0; i < a; ++i) {
        std::cout << "----------";
    }
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << std::fixed <<std::setprecision(4) <<std::setw(7)<< " || "<<
              std::setw(6) << 'h' <<" | " <<
              std::setw(7) << "dis" << " | ";

    for(size_t i = 0; i < a; ++i){
        if(i==a/2) {
            std::cout << std::fixed << std::setprecision(4) <<
                      std::setw(7) << "  alpha ";
        }
        else std::cout << std::fixed << std::setprecision(4) <<
                       std::setw(7) << "           ";
    }

    std::cout << std::fixed <<std::setprecision(4) <<
              std::setw(7) << " |     w  " << " | " <<
              std::setw(7) << "     d" << " || "<< std::endl<<"    ";
    for(size_t i = 0; i < a; ++i) {
        std::cout << "----------";
    }
    std::cout << "-----------------------------------------------" << std::endl;
}

void disp2(double h, double dis, std::vector <double> alpha, double w, double d){

    std::cout << std::fixed << std::setprecision(4) <<
                       std::setw(7) << " || " <<
                       std::setw(2) << h << " | " <<
                       std::setw(7) << dis << " | [  ";

    for(size_t i = 0; i < alpha.size(); ++i){
        std::cout << std::fixed << std::setprecision(4) <<
                  std::setw(7) << alpha[i] << "  ";
    }

    std::cout   <<
                       std::setw(7) << "] | " << std::fixed << std::setprecision(4)<<
                       std::setw(7) << w << " | " <<
                       std::setw(7) << d << " ||"<< std::endl << "    -";
    for(size_t i = 0; i < alpha.size(); ++i) {
        std::cout << "--------";
    }
    std::cout << "------------------------------------------------------" << std::endl;

}



double srGeometr(std::vector <double> fun1, std::vector <double> fun2){
    double sum = 1;
    for(int i = 0; i < fun2.size(); ++i){
        sum*=pow(fun1[i], (fun2[i]));
    }
    return sum;//fun1.size();
}

double cheb1(std::vector <double> fun){
    double max = fun[1]-fun[0];
    for(int i = 2; i < 100; ++i)
        if((fun[i]-fun[i-1])*(fun[i]-fun[i-1]) > max*max)
            max = pow(((fun[i]-fun[i-1])*(fun[i]-fun[i-1])),(1./2));
    return max;
}

double cheb2(std::vector <double> fun1, std::vector <double> fun2){
    double max = fun1[0]-fun2[0];
    for(int i = 1; i < 100; ++i){
        if((fun1[i]-fun2[i])*(fun1[i]-fun2[i]) > max*max){
            max = pow(((fun1[i]-fun2[i])*(fun1[i]-fun2[i])),(1./2));
        }
    }
    return max;
}

int main() {
    std::cout << "Hello, World!" << std::endl;

    double xMin = 0;
    double xMax = 3.14;
    double ample = 0.25;
    double P = 0.95;
    double eps = 0.01;

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution <double> dis(-ample,ample);

    double K = 100; // количество отсчетов

    std::vector <double> xk; // наши иксы
    for(int i = 0; i <= K; ++i){
        double one = xMin + i*(xMax - xMin)/K;
        xk.push_back(one);
    }
    std::vector <double> fxk; // значения функции
    std::vector <double> fxk_noised; // значение зашумленной функции

    for(int i = 0; i <= K; ++i){
        double one = function(xk[i]);
        fxk.push_back(one);
        fxk_noised.push_back(one+dis(gen));
    }


    std::vector <double> r { 3.0 , 5.0 }; // размер скользящего окна
    double r_f = 3;
    double r_s = 5;

    double L = 10;
    std::vector <double> h;

    for(double i = 0.0; i <= 1; i+=0.1){
        h.push_back(i);
    }
    std::vector <double> lambda; // значения весов
    for(size_t i=0;i<=10;++i){
        double one = i / L;
        lambda.push_back(one);
    }


    std::uniform_real_distribution<double> dis2(0.000,1.000);
    double N = log(1. - P)/log(1. - eps/(xMax - xMin));

    int m_1 = (r[0]+1.0)/2.0;
    int m_2 = (r[1]+1.0)/2.0;
    Answer best;
    std::vector <Answer> answers1;

    std::cout << "  r = 3 "<<std::endl;
    head(3);
//3
    for(size_t i = 0 ; i < h.size();++i) {
        Answer best_p;

        for(size_t j = 0; j < N; ++j) {

            std::vector <double> alpha;
            alpha.push_back(0);
            alpha.push_back(0);
            alpha.push_back(0);

            alpha[1] = std::uniform_real_distribution<double>{0, 1}(gen);
            alpha[0] = alpha[2] = 0.5 * (1 - alpha[1]);


            std::vector <double> fxk_filtered; // Вектор значений отфильтрованной функции
            for(size_t t = 0; t < K; ++t){
                if(t==0 || t== K-1) {
                    fxk_filtered.push_back(fxk_noised[t]);
                }
                else {
                    std::vector <double> fxk_npie;
                    fxk_npie.push_back(fxk_noised[t-1]);
                    fxk_npie.push_back(fxk_noised[t]);
                    fxk_npie.push_back(fxk_noised[t+1]);
                    double ons = srGeometr(fxk_npie, alpha);
                    fxk_filtered.push_back(ons);
                }
            }


            double w = cheb1(fxk_filtered);
            double d = cheb2(fxk_filtered, fxk_noised);

            double J = w*lambda[i]+(1-lambda[i])*d;

            double dis = w > d ? w : d;

            if(j == 0 || J > best.dis){
                best_p.lambda=lambda[i];//==h
                best_p.dis=dis;//==final criterij
                best_p.alpha=alpha;
                best_p.w=w;
                best_p.d=d;
                best_p.J=J;
            }
        }


        disp2(best_p.lambda, best_p.dis, best_p.alpha, best_p.w, best_p.d);

        if(i==0 || best.dis > best_p.dis){
            best = best_p;
        }
    }

    disp1(best.lambda, best.J, best.w, best.d);

//======================================================================================================================
    for(size_t i = 0; i < K; ++i){
        std::cout<<'('<< xk[i] <<';'<< fxk[i] <<')';
    }
    std::cout<< std::endl;

    for(size_t i = 0; i < K; ++i){
        std::cout<<'('<< xk[i] <<';'<< fxk_noised[i] <<')';
    }
    std::cout<< std::endl;
    std::vector <double> fxk_filtered1; // Вектор значений отфильтрованной функции
    for(size_t t = 0; t < K; ++t){
        if(t==0 || t== K-1) {
            fxk_filtered1.push_back(fxk_noised[t]);
        }
        else {
            std::vector <double> fxk_npie1;
            fxk_npie1.push_back(fxk_noised[t-1]);
            fxk_npie1.push_back(fxk_noised[t]);
            fxk_npie1.push_back(fxk_noised[t+1]);
            double ons = srGeometr(fxk_npie1, best.alpha);
            fxk_filtered1.push_back(ons);
        }
    }
    for(size_t i = 0; i < K; ++i){
        std::cout <<'('<< xk[i] <<';'<< fxk_filtered1[i] <<')';
    }
    std::cout<< std::endl;
//======================================================================================================================

    std::vector <Answer> answers2;

    std::cout << "  r = 5 "<<std::endl;
    head2(5);
//5
    for(size_t i = 0 ; i < h.size();++i) {
        Answer best_p;

        for (size_t j = 0; j < N; ++j) {

            std::vector<double> alpha;
            alpha.push_back(0);
            alpha.push_back(0);
            alpha.push_back(0);
            alpha.push_back(0);
            alpha.push_back(0);

            alpha[2] = std::uniform_real_distribution<double>{0, 1}(gen);
            alpha[1] = alpha[3] = 0.3 * (1 - alpha[2]);
            alpha[0] = alpha[4] = 0.5 * (alpha[1]);

            std::vector<double> fxk_filtered; // Вектор значений отфильтрованной функции
            for (size_t t = 0; t < K; ++t) {
                if (t == 0 || t==1 || t == K - 1 || t==K-2) {
                    fxk_filtered.push_back(fxk_noised[t]);
                } else {
                    std::vector<double> fxk_npie;
                    fxk_npie.push_back(fxk_noised[t - 2]);
                    fxk_npie.push_back(fxk_noised[t - 1]);
                    fxk_npie.push_back(fxk_noised[t]);
                    fxk_npie.push_back(fxk_noised[t + 1]);
                    fxk_npie.push_back(fxk_noised[t + 2]);
                    double ons = srGeometr(fxk_npie, alpha);
                    fxk_filtered.push_back(ons);
                }
            }


            double w = cheb1(fxk_filtered);
            double d = cheb2(fxk_filtered, fxk_noised);

            double J = w * lambda[i] + (1 - lambda[i]) * d;

            double dis = w > d ? w : d;

            if (j == 0 || J > best.dis) {
                best_p.lambda = lambda[i];//==h
                best_p.dis = dis;//==final criterij
                best_p.alpha = alpha;
                best_p.w = w;
                best_p.d = d;
                best_p.J = J;
            }
        }


        disp2(best_p.lambda, best_p.dis, best_p.alpha, best_p.w, best_p.d);

        if (i == 0 || best.dis > best_p.dis) {
            best = best_p;
        }
    }

    disp1(best.lambda, best.J, best.w, best.d);
    // вывод



    std::cout<< std::endl;
    std::vector<double> fxk_filtered2; // Вектор значений отфильтрованной функции
    for (size_t t = 0; t < K; ++t) {
        if (t == 0 || t==1 || t == K - 1 || t==K-2) {
            fxk_filtered2.push_back(fxk_noised[t]);
        } else {
            std::vector<double> fxk_npie2;
            fxk_npie2.push_back(fxk_noised[t - 2]);
            fxk_npie2.push_back(fxk_noised[t - 1]);
            fxk_npie2.push_back(fxk_noised[t]);
            fxk_npie2.push_back(fxk_noised[t + 1]);
            fxk_npie2.push_back(fxk_noised[t + 2]);
            double ons = srGeometr(fxk_npie2, best.alpha);
            fxk_filtered2.push_back(ons);
        }
    }
    for(size_t i = 0; i < K; ++i){
        std::cout <<'('<< xk[i] <<';'<< fxk_filtered2[i] <<')';
    }
    return 0;
}