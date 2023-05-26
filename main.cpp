// Dimensiunea populatiei:
// Domeniul de definitie al functiei(interval inchis):
// Parametrii pentru functia de maximizat(coeficientii):
// Precizia cu care se lucreaza:
// Probabilitatea de recombinare:
// Probabilitatea de mutatie:
// Numarul de etape al algoritmului:

#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;

ifstream fin("input.in");
ofstream fout("output.out");

int cautareBin(vector<double>& intervale, double u) {
    int st = 0;
    int dr = intervale.size() - 1;
    while (st <= dr) {
        int m = st + (dr - st) / 2;
        if (intervale[m] <= u && (m == intervale.size() - 1 || u < intervale[m + 1])) {
            return m;
        }
        else if (intervale[m] < u) {
            st = m + 1;
        }
        else {
            dr = m - 1;
        }
    }
    return -1;
}


int main() {
    /// initializare si citire
    int dim, etape;
    double prec, a, b, x, y, z, probRecomb, probMut;
    srand(time(NULL));
    fin>>dim>>a>>b>>x>>y>>z>>prec>>probRecomb>>probMut>>etape;
    vector<vector<int>> binPop(dim);
    vector<double> valx(dim), valf(dim);
    int random;
    int l = ceil(log2(double((b - a) * pow(10, prec))));
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < l; ++j) {
            binPop[i].push_back(rand() % 2);
        }
    }
    for (int i = 0; i < dim; ++i) {
        double aux = 0, power = l - 1;
        for (int j = 0; j < binPop[i].size(); ++j) {
            aux = aux + binPop[i][j] * pow(2, power);
            power--;
        }
        valx[i] =((b - a) / (pow(2, l) - 1)) * aux + a;
        valf[i] = valx[i] * valx[i] * x + valx[i] * y + z;
    }

    double bestFit1 = INT_MIN, bestFitx;
    vector<int> bestFit;
    for (int i = 0; i < dim; ++i) {
        if (bestFit1 < valf[i]){
            bestFit1 = valf[i];
            bestFitx = valx[i];
            bestFit = binPop[i];
        }
    }

    fout<<"Populatia initiala"<<endl;
    for (int i = 0; i < binPop.size(); ++i) {
        fout<<i + 1<<": ";
        for (int j = 0; j < binPop[i].size(); ++j) {
            fout<<binPop[i][j];
        }
        fout<<"   x: "<<fixed<<setprecision(prec)<<valx[i];
        fout<<"   f: "<<fixed<<setprecision(15)<<valf[i];
        fout<<endl;
    }
    fout<<endl;
    double sumf = 0;
    for (int i = 0; i < dim; ++i) {
        sumf = sumf + valf[i];
    }
    vector<double> intervale(dim + 1);
    intervale[0] = 0;



    fout<<"Probabilitati selectie"<<endl;
    for (int i = 0; i < dim; ++i) {
        intervale[i + 1] = valf[i] / sumf + intervale[i];
        fout<<"cromozom "<<i + 1<<" probabilitate "<<fixed<<setprecision(15)<<valf[i] / sumf<<endl;
    }
    fout<<endl;


    fout<<"Intervale probabilitati selectie"<<endl;
    for (int i = 0; i < dim; ++i) {
        fout<<i + 1<<": "<<fixed<<setprecision(15)<<intervale[i]<<' '<<intervale[i + 1]<<endl;
    }
    fout<<endl;


    vector<vector<int>> aux(dim);
    int k = -1;
    for (int i = 0; i < dim; ++i) {
        double u = 0;
        for (int j = 0; j < 14; ++j) {
            int aux = rand();
            u = u/10 + aux%10;
        }
        u = u / 10;
        int index = cautareBin(intervale, u);
        aux[++k] = binPop[index];
        fout<<"u = "<<fixed<<setprecision(15)<<u<<" selectam cromozomul "<<index + 1<<" (noul cromozom "<<k + 1<<')';
        fout<<endl;
    }
    fout<<endl;


    fout<<"Dupa selectie"<<endl;
    binPop = aux;
    for (int i = 0; i < dim; ++i) {
        double aux = 0, power = l - 1;
        for (int j = 0; j < binPop[i].size(); ++j) {
            aux = aux + binPop[i][j] * pow(2, power);
            power--;
        }
        valx[i] = ((b - a) / (pow(2, l) - 1)) * aux + a;
        valf[i] = valx[i] * valx[i] * x + valx[i] * y + z;
    }
    for (int i = 0; i < dim; ++i) {
        fout<<i + 1<<": ";
        for (int j = 0; j < binPop[i].size(); ++j) {
            fout<<binPop[i][j];
        }
        fout<<"   x: "<<fixed<<setprecision(prec)<<valx[i];
        fout<<"   f: "<<fixed<<setprecision(15)<<valf[i];
        fout<<endl;
    }
    fout<<endl;


    vector<pair<double, int>> participa;
    fout<<"Probabilitatea de incrucisare "<<probRecomb<<endl;
    for (int i = 0; i < dim; ++i) {
        fout<<i + 1<<": ";
        for (int j = 0; j < binPop[i].size(); ++j) {
            fout<<binPop[i][j];
        }
        double u = 0;
        for (int j = 0; j < 14; ++j) {
            int aux = rand();
            u = u/10 + aux%10;
        }
        u = u / 10;
        fout<<" u = "<<u;
        if (u < probRecomb){
            fout<<" < "<<probRecomb<<" participa";
            participa.push_back(make_pair(valx[i], i));
        }
        fout<<endl;
    }
    sort(participa.begin(), participa.end());
    reverse(participa.begin(), participa.end());


    fout<<endl;
    for (int i = 0; i < participa.size() - 1; ++i) {
        if (participa[i].second == -1)
            continue;
        int j = i + 1;
        while(participa[j].first == participa[i].first
              && j < participa.size()
              && participa[j].second == -1){
            ++j;
        }
        if (j == participa.size() - 1){
            if (participa[i].first != participa[j].first && participa[j].second != -1){
                fout<<"Recombinare dintre cromozomul "<<participa[i].second + 1<<" cu cromozomul "<<participa[j].second + 1<<":\n";
                for (int ll = 0; ll < binPop[participa[i].second].size(); ++ll) {
                    fout<<binPop[participa[i].second][ll];
                }
                fout<<' ';
                for (int ll = 0; ll < binPop[participa[j].second].size(); ++ll) {
                    fout<<binPop[participa[j].second][ll];
                }
                int punct = rand() % l;
                fout<<" punct "<<punct<<endl;
                fout<<"Rezultat ";
                for (int ll = punct; ll < l; ++ll) {
                    int aux = binPop[participa[j].second][ll];
                    binPop[participa[j].second][ll] = binPop[participa[i].second][ll];
                    binPop[participa[i].second][ll] = aux;
                }
                for (int ll = 0; ll < binPop[participa[i].second].size(); ++ll) {
                    fout<<binPop[participa[i].second][ll];
                }
                fout<<' ';
                for (int ll = 0; ll < binPop[participa[j].second].size(); ++ll) {
                    fout<<binPop[participa[j].second][ll];
                }
                fout<<endl<<endl;
            }
            break;
        }
        else{
            fout<<"Recombinare dintre cromozomul "<<participa[i].second + 1<<" cu cromozomul "<<participa[j].second + 1<<":\n";
            for (int ll = 0; ll < binPop[participa[i].second].size(); ++ll) {
                fout<<binPop[participa[i].second][ll];
            }
            fout<<' ';
            for (int ll = 0; ll < binPop[participa[j].second].size(); ++ll) {
                fout<<binPop[participa[j].second][ll];
            }
            int punct = rand() % l;
            fout<<" punct "<<punct<<endl;
            fout<<"Rezultat ";
            for (int ll = punct; ll < l; ++ll) {
                int aux = binPop[participa[j].second][ll];
                binPop[participa[j].second][ll] = binPop[participa[i].second][ll];
                binPop[participa[i].second][ll] = aux;
            }
            for (int ll = 0; ll < binPop[participa[i].second].size(); ++ll) {
                fout<<binPop[participa[i].second][ll];
            }
            fout<<' ';
            for (int ll = 0; ll < binPop[participa[j].second].size(); ++ll) {
                fout<<binPop[participa[j].second][ll];
            }
            double aux = 0, power = l - 1;
            for (int ll = 0; ll < binPop[participa[i].second].size(); ++ll) {
                aux = aux + binPop[participa[i].second][ll] * pow(2, power);
                power--;
            }
            valx[participa[i].second] =((b - a) / (pow(2, l) - 1)) * aux + a;
            valf[participa[i].second] = valx[participa[i].second] * valx[participa[i].second] * x + valx[participa[i].second] * y + z;
            aux = 0;
            power = l - 1;
            for (int ll = 0; ll < binPop[participa[j].second].size(); ++ll) {
                aux = aux + binPop[participa[j].second][ll] * pow(2, power);
                power--;
            }
            valx[participa[j].second] =((b - a) / (pow(2, l) - 1)) * aux + a;
            valf[participa[j].second] = valx[participa[j].second] * valx[participa[j].second] * x + valx[participa[j].second] * y + z;
            participa[i].second = -1;
            participa[j].second = -1;
            fout<<endl<<endl;
        }
    }
    fout<<endl;


    fout<<"Dupa recombinare"<<endl;
    for (int i = 0; i < binPop.size(); ++i) {
        fout<<i + 1<<": ";
        for (int j = 0; j < binPop[i].size(); ++j) {
            fout<<binPop[i][j];
        }
        fout<<"   x: "<<fixed<<setprecision(prec)<<valx[i];
        fout<<"   f: "<<fixed<<setprecision(15)<<valf[i];
        fout<<endl;
    }
    fout<<endl;


    fout<<"Probabilitate de mutatie pentru fiecare gena "<<probMut<<endl;
    vector<int> partMut;
    int ok = 0;
    for (int i = 0; i < dim; ++i) {
        int ok1 = 0;
        for (int ll = 0; ll < binPop[i].size(); ++ll) {
            double u = 0;
            for (int j = 0; j < 14; ++j) {
                int aux = rand();
                u = u/10 + aux%10;
            }
            u = u / 10;
            if (u < probMut){
                ok = ok1 = 1;
                binPop[i][ll] = 1 - binPop[i][ll];
            }
        }
        if (ok1 == 1){
            partMut.push_back(i);
            double aux = 0, power = l - 1;
            for (int j = 0; j < binPop[i].size(); ++j) {
                aux = aux + binPop[i][j] * pow(2, power);
                power--;
            }
            valx[i] =((b - a) / (pow(2, l) - 1)) * aux + a;
            valf[i] = valx[i] * valx[i] * x + valx[i] * y + z;
        }
    }

    if (ok == 1){
        fout<<"Au fost modificati cromozomii: "<<endl;
        for (int i = 0; i < partMut.size(); ++i) {
            fout<<partMut[i] + 1<<endl;
        }
    }
    fout<<endl;


    fout<<"Dupa mutatie"<<endl;
    for (int i = 0; i < binPop.size(); ++i) {
        fout<<i + 1<<": ";
        for (int j = 0; j < binPop[i].size(); ++j) {
            fout<<binPop[i][j];
        }
        fout<<"   x: "<<fixed<<setprecision(prec)<<valx[i];
        fout<<"   f: "<<fixed<<setprecision(15)<<valf[i];
        fout<<endl;
    }
    fout<<endl;

    for (int i = 0; i < dim; ++i) {
        if (bestFit1 < valf[i]){
            bestFit1 = valf[i];
            bestFitx = valx[i];
            bestFit = binPop[i];
        }
    }


    fout<<"Evolutia maximului si evolutia mediei"<<endl;
    double maxf = INT_MIN, medf = 0;
    for (int i = 0; i < dim; ++i) {
        if (maxf < valf[i]){
            maxf = valf[i];
        }
        medf = medf + valf[i];
    }
    medf = medf / dim;
    fout<<fixed<<setprecision(15)<<maxf<<"    "<<medf<<endl;
    for (int i = 0; i < etape - 1; ++i) {
        for (int i = 0; i < dim; ++i) {
            if (bestFit1 < valf[i]){
                bestFit1 = valf[i];
                bestFitx = valx[i];
                bestFit = binPop[i];
            }
        }
        for (int i = 0; i < dim; ++i) {
            double aux = 0, power = l - 1;
            for (int j = 0; j < binPop[i].size(); ++j) {
                aux = aux + binPop[i][j] * pow(2, power);
                power--;
            }
            valx[i] = ((b - a) / (pow(2, l) - 1)) * aux + a;
            valf[i] = valx[i] * valx[i] * x + valx[i] * y + z;
        }


        double sumf = 0;
        for (int i = 0; i < dim; ++i) {
            sumf = sumf + valf[i];
        }
        vector<double> intervale(dim + 1);
        intervale[0] = 0;


        for (int i = 0; i < dim; ++i) {
            intervale[i + 1] = valf[i] / sumf + intervale[i];
        }


        vector<vector<int>> aux(dim);
        int k = -1;
        for (int i = 0; i < dim; ++i) {
            double u = 0;
            for (int j = 0; j < 14; ++j) {
                int auxx = rand();
                u = u/10 + auxx%10;
            }
            u = u / 10;
            int index = cautareBin(intervale, u);
            aux[++k] = binPop[index];
        }
        binPop = aux;
        for (int i = 0; i < dim; ++i) {
            double auxx = 0, power = l - 1;
            for (int j = 0; j < binPop[i].size(); ++j) {
                auxx = auxx + binPop[i][j] * pow(2, power);
                power--;
            }
            valx[i] = ((b - a) / (pow(2, l) - 1)) * auxx + a;
            valf[i] = valx[i] * valx[i] * x + valx[i] * y + z;
        }
        vector<pair<double, int>> participa;
        for (int i = 0; i < dim; ++i) {
            double u = 0;
            for (int j = 0; j < 14; ++j) {
                int aux = rand();
                u = u/10 + aux%10;
            }
            u = u / 10;
            if (u < probRecomb){
                participa.push_back(make_pair(valx[i], i));
            }
        }
        sort(participa.begin(), participa.end());
        reverse(participa.begin(), participa.end());
        for (int i = 0; i < participa.size() - 1; ++i) {
            if (participa[i].second == -1)
                continue;
            int j = i + 1;
            while(participa[j].first == participa[i].first
                  && j < participa.size()
                  && participa[j].second == -1){
                ++j;
            }
            if (j == participa.size() - 1){
                if (participa[i].first != participa[j].first && participa[j].second != -1){
                    int punct = rand() % l;
                    for (int ll = punct; ll < l; ++ll) {
                        int aux = binPop[participa[j].second][ll];
                        binPop[participa[j].second][ll] = binPop[participa[i].second][ll];
                        binPop[participa[i].second][ll] = aux;
                    }
                }
                break;
            }
            else{
                int punct = rand() % l;
                for (int ll = punct; ll < l; ++ll) {
                    int aux = binPop[participa[j].second][l];
                    binPop[participa[j].second][ll] = binPop[participa[i].second][ll];
                    binPop[participa[i].second][ll] = aux;
                }
                double aux = 0, power = l - 1;
                for (int ll = 0; ll < binPop[participa[i].second].size(); ++ll) {
                    aux = aux + binPop[participa[i].second][ll] * pow(2, power);
                    power--;
                }
                valx[participa[i].second] = ((b - a) / (pow(2, l) - 1)) * aux + a;
                valf[participa[i].second] = valx[participa[i].second] * valx[participa[i].second] * x + valx[participa[i].second] * y + z;
                aux = 0;
                power = l - 1;
                for (int ll = 0; ll < binPop[participa[j].second].size(); ++ll) {
                    aux = aux + binPop[participa[j].second][ll] * pow(2, power);
                    power--;
                }
                valx[participa[j].second] = ((b - a) / (pow(2, l) - 1)) * aux + a;
                valf[participa[j].second] = valx[participa[j].second] * valx[participa[j].second] * x + valx[participa[j].second] * y + z;
                participa[i].second = -1;
                participa[j].second = -1;
            }
        }
        vector<int> partMut;
        for (int i = 0; i < dim; ++i) {
            int ok1 = 0;
            for (int ll = 0; ll < binPop[i].size(); ++ll) {
                double u = 0;
                for (int j = 0; j < 14; ++j) {
                    int aux = rand();
                    u = u/10 + aux%10;
                }
                u = u / 10;
                if (u < probMut){
                    ok1 = 1;
                    binPop[i][ll] = 1 - binPop[i][ll];
                }
            }
            if (ok1 == 1){
                double aux = 0, power = l - 1;
                for (int j = 0; j < binPop[i].size(); ++j) {
                    aux = aux + binPop[i][j] * pow(2, power);
                    power--;
                }
                valx[i] = ((b - a) / (pow(2, l) - 1)) * aux + a;
                valf[i] = valx[i] * valx[i] * x + valx[i] * y + z;
            }
        }

        int okbest = 0;
        for (int j = 0; j < dim; ++j) {
            if (bestFit1 < valf[i]){
                bestFit1 = valf[i];
                bestFitx = valx[i];
                bestFit = binPop[i];
                okbest  = 1;
            }
            else if (bestFit1 == valf[i]){
                okbest = 1;
            }
        }

        double worst = INT_MAX;
        int worstp;
        if (okbest == 0){
            for (int j = 0; j < dim; ++j) {
                if (worst > valf[i]){
                    worst = valf[i];
                    worstp = j;
                }
            }
            valf[worstp] = bestFit1;
            valx[worstp] = bestFitx;
            binPop[worstp] = bestFit;
        }

        for (int i = 0; i < dim; ++i) {
            if (maxf < valf[i]){
                maxf = valf[i];
            }
            medf = medf + valf[i];
        }
        medf = medf / dim;
        fout<<fixed<<setprecision(15)<<maxf<<"    "<<medf<<endl;
    }
    return 0;
}
