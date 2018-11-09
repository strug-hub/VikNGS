#include "Math.h"
#include <random>

static std::random_device rd;
static thread_local std::mt19937 generate(rd());

int randomInt(int from, int to) {
    std::uniform_int_distribution<> sample(from, to);
    return sample(generate);
}
double randomDouble(double from, double to) {
    std::uniform_real_distribution<double> sample(from, to);
    return sample(generate);
}
double randomNormal(double mean, double sd) {
    std::normal_distribution<> sample(mean, sd);
    return sample(generate);
}
int randomBinomial(int trials, double success) {
    std::binomial_distribution<> sample(trials, success);
    return sample(generate);
}

MatrixXd shuffleColumnwiseWithoutReplacement(MatrixXd& M){

    MatrixXd shuffled(M.rows(), M.cols());
    VectorXi indices = VectorXi::LinSpaced(M.rows(), 0, M.rows());
    for(int i = 0; i < M.cols(); i++){
        std::random_shuffle(indices.data(), indices.data() + M.rows());
        shuffled.col(i) = indices.asPermutation() * M.col(i);
    }

    return shuffled;
}

VectorXd shuffleWithoutReplacement(VectorXd& V){

    VectorXd shuffled(V.rows());
    VectorXi indices = VectorXi::LinSpaced(V.rows(), 0, V.rows());
    std::random_shuffle(indices.data(), indices.data() + V.rows());
    shuffled = indices.asPermutation() * V;

    return shuffled;
}

MatrixXd groupwiseShuffleWithReplacement(MatrixXd& M, VectorXi& G, std::map<int, std::vector<int>>& group){

    MatrixXd shuffled(M.rows(), M.cols());
    int g, rand;

    for(int i = 0; i < M.rows(); i++){

        g = G[i];
        rand = randomInt(0, group[g].size() - 1);
        rand = group[g][static_cast<size_t>(rand)];

        shuffled.row(i) = M.row(rand);
    }

    return shuffled;
}

VectorXd groupwiseShuffleWithReplacement(VectorXd& V, VectorXi& G, std::map<int, std::vector<int>>& group){

    VectorXd shuffled(V.rows());
    int g, n, rand;

    for(int j = 0; j < V.rows(); j++){
        g = G[j]; n = static_cast<int>(group[g].size());
        rand = randomInt(0, n - 1);
        rand = group[g][static_cast<size_t>(rand)];

        shuffled(j) = V[rand];
    }

    return shuffled;
}

VectorXd groupwiseShuffleWithoutReplacement(VectorXd& V, VectorXi& G, std::map<int, std::vector<int>>& group){

    std::map<int, std::vector<int>> groupShuffle;

    for (std::pair<int, std::vector<int>> e : group){
        std::vector<int> v = e.second;
        std::random_shuffle(v.begin(), v.end());
        groupShuffle[e.first] = v;
    }

    VectorXd shuffled(V.rows());
    int g;

    for(int j = 0; j < V.rows(); j++){
        g = G[j];

        shuffled(j) = V[groupShuffle[g].back()];
        groupShuffle[g].pop_back();
    }

    return shuffled;
}



