#include "Math.h"
#include <random>

static std::random_device rd;
static std::mt19937 generate(rd());

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

std::vector<VectorXd> groupwiseShuffleWithoutReplacement(std::vector<VectorXd> &v){

    std::vector<VectorXd> shuffled;

    for(size_t i = 0; i < v.size(); i++){

        VectorXi indices = VectorXi::LinSpaced(v[i].rows(), 0, v[i].rows());
        std::random_shuffle(indices.data(), indices.data() + v[i].rows());
        VectorXd s = indices.asPermutation() * v[i];

        shuffled.push_back(s);
    }

    return shuffled;
}

std::vector<VectorXd> groupwiseShuffleWithReplacement(std::vector<VectorXd> &v){

    std::vector<VectorXd> shuffled;
    int n, j;

    for(size_t i = 0; i < v.size(); i++){
        n = v[i].rows();
        VectorXd s(n);

        for (j = 0; j < n; j++)
            s[j] = v[i][randomInt(0, n - 1)];

        shuffled.push_back(s);
    }

    return shuffled;

}

std::vector<MatrixXd> groupwiseShuffleWithReplacement(std::vector<MatrixXd> &m){

    std::vector<MatrixXd> shuffled;
    int n, ncol, j, k;

    for(size_t i = 0; i < m.size(); i++){
        n = m[i].rows();
        ncol = m[i].cols();
        MatrixXd s(n, ncol);

        for (k = 0; k < ncol; k++)
            for (j = 0; j < n; j++)
                s(j,k) = m[i](randomInt(0, n - 1),k);

        shuffled.push_back(s);
    }

    return shuffled;

}

std::vector<MatrixXd> shuffleWithoutReplacement(std::vector<MatrixXd> &m){

    MatrixXd cat = concatenate(m);

    VectorXi indices = VectorXi::LinSpaced(cat.rows(), 0, cat.rows());
    std::random_shuffle(indices.data(), indices.data() + cat.rows());
    MatrixXd s = indices.asPermutation() * cat;

    std::vector<MatrixXd> shuffled;
    int startRow = 0;

    for(size_t i = 0; i < m.size(); i++){
        MatrixXd m_i = s.block(0, startRow, m[i].rows(), m[i].cols());
        startRow += m[i].rows();
        shuffled.push_back(m_i);
    }

    return shuffled;

}

std::vector<VectorXd> shuffleWithoutReplacement(std::vector<VectorXd> &v){

    VectorXd cat = concatenate(v);

    VectorXi indices = VectorXi::LinSpaced(cat.rows(), 0, cat.rows());
    std::random_shuffle(indices.data(), indices.data() + cat.rows());
    VectorXd s = indices.asPermutation() * cat;

    std::vector<VectorXd> shuffled;
    int startRow = 0;

    for(size_t i = 0; i < v.size(); i++){
        VectorXd v_i = s.segment(startRow, v[i].rows());
        startRow += v[i].rows();
        shuffled.push_back(v_i);
    }

    return shuffled;

}
