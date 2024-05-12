#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const long double tol = 1e-8;
int maxIter = 2000;
 
long double emission_prob(long double mean, long double variance, long double y_t) {
    const long double pi = M_PI;
    long double r = (1 / sqrt(2 * pi * variance)) * (exp(-pow(y_t - mean, 2) / (2 * variance)));
    return r;
}

vector<vector<long double>> forward_prob(vector<vector<long double>>& a, vector<long double>& initial_prob, int T, vector<long double>& observation, vector<long double>& mean, vector<long double>& variance) {
    vector<vector<long double>> res(T, vector<long double>(2));
    res[0][0] = initial_prob[0] * emission_prob(mean[0], variance[0], observation[0]);
    res[0][1] = initial_prob[1] * emission_prob(mean[1], variance[1], observation[0]);
    for (int i = 1; i < T; i++) {
        res[i][0] = emission_prob(mean[0], variance[0], observation[i]) * ((a[0][0] * res[i - 1][0]) + (a[1][0] * res[i - 1][1]));
        res[i][1] = emission_prob(mean[1], variance[1], observation[i]) * ((a[0][1] * res[i - 1][0]) + (a[1][1] * res[i - 1][1]));
    }
    return res;
}

vector<vector<long double>> backward_prob(vector<vector<long double>>& a, vector<long double>& initial_prob, int T, vector<long double>& observation, vector<long double>& mean, vector<long double>& variance) {
    vector<vector<long double>> res(T, vector<long double>(2));
    res[T - 1][0] = res[T - 1][1] = 1; 
    for (int i = T - 2; i >= 0; i--) {
        res[i][0] = (a[0][0] * emission_prob(mean[0], variance[0], observation[i + 1]) * res[i + 1][0]) + (a[0][1] * emission_prob(mean[1], variance[1], observation[i + 1]) * res[i + 1][1]);
        res[i][1] = (a[1][0] * emission_prob(mean[0], variance[0], observation[i + 1]) * res[i + 1][0]) + (a[1][1] * emission_prob(mean[1], variance[1], observation[i + 1]) * res[i + 1][1]);
    }
    return res;
}

vector<vector<long double>> state_posterior_prob(vector<vector<long double>>& f, int T, vector<vector<long double>>& b) {
    vector<vector<long double>> res(T, vector<long double>(2));
    for (int i = 0; i < T; i++) {
        long double denominator = (f[i][0] * b[i][0]) + (f[i][1] * b[i][1]);
        res[i][0] = (f[i][0] * b[i][0]) / denominator;
        res[i][1] = (f[i][1] * b[i][1]) / denominator;
    }
    return res;
}

vector<vector<vector<long double>>> transition_posterior_prob(vector<vector<long double>>& alpha, vector<vector<long double>>& a, int T, vector<vector<long double>>& beta, vector<long double>& mean, vector<long double>& variance, vector<long double>& observation) {
    vector<vector<vector<long double>>> epsilon(T - 1, vector<vector<long double>>(2, vector<long double>(2)));

    for (int t = 0; t < T - 1; ++t) {
        long double den = 0.0;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                den += alpha[t][i] * a[i][j] * emission_prob(mean[j], variance[j], observation[t + 1]) * beta[t + 1][j];
            }
        }
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                long double num = alpha[t][i] * a[i][j] * emission_prob(mean[j], variance[j], observation[t + 1]) * beta[t + 1][j];
                epsilon[t][i][j] = num / den;
            }
        }
    }
    return epsilon;
}

int main() {
    vector<vector<long double>> trans_matrix(2, vector<long double>(2));
    long double t1, t2, t3, t4;
    cin >> t1 >> t2 >> t3 >> t4;
    trans_matrix[0][0] = t1;
    trans_matrix[0][1] = t2;
    trans_matrix[1][0] = t3;
    trans_matrix[1][1] = t4;

    vector<long double> initial_prob(2);
    long double p1, p2;
    cin >> p1 >> p2;
    initial_prob[0] = p1;
    initial_prob[1] = p2;

    int T;
    cin >> T;

    vector<long double> observation(T);
    for (int i = 0; i < T; i++) {
        cin >> observation[i];
    }
    // auto start = std::chrono::high_resolution_clock::now();
    
    vector<long double> mean = {0,0};
    vector<long double> variance = {1,1};

    vector<long double> new_mean(2);
    vector<long double> new_variance(2);
    vector<long double> new_initial_prob(2, 0.5); 
    vector<vector<long double>> new_trans_matrix(2, vector<long double>(2, 0.5)); // Initialize with 0.5 for each transition

    bool converged = false;
    int iter = 0;
    
    while (!converged && iter<maxIter) {
        vector<vector<long double>> alpha = forward_prob(trans_matrix, initial_prob, T, observation, mean, variance);
        vector<vector<long double>> beta = backward_prob(trans_matrix, initial_prob, T, observation, mean, variance);
        vector<vector<long double>> gamma = state_posterior_prob(alpha, T, beta);
        vector<vector<vector<long double>>> epsilon = transition_posterior_prob(alpha, trans_matrix, T, beta, mean, variance, observation);

        // Update initial probabilities
        new_initial_prob[0] = gamma[0][0];
        new_initial_prob[1] = gamma[0][1];

        // Update transition matrix
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                long double num = 0;
                long double den = 0;
                for (int t = 0; t < T - 1; t++) {
                    num += epsilon[t][i][j];
                    den += gamma[t][i]; 
                }
                new_trans_matrix[i][j] = num / den;
            }
        }

        // Update mean and variance
        long double sum1 = 0;
        long double sum2 = 0;
        long double var1 = 0, var2 = 0;
        long double den_m1 = 0;
        long double den_m2 = 0;
        for (int i = 0; i < T; i++) {
            sum1 += observation[i] * gamma[i][0];
            sum2 += observation[i] * gamma[i][1];
            var1 += gamma[i][0] * pow(observation[i] - mean[0], 2);
            var2 += gamma[i][1] * pow(observation[i] - mean[1], 2);
            den_m1 += gamma[i][0];
            den_m2 += gamma[i][1];
        }
        new_mean[0] = sum1 / den_m1;
        new_mean[1] = sum2 / den_m2;

        new_variance[0] = var1/den_m1;
        new_variance[1] = var2/den_m2;

        long double piDiff = abs(new_initial_prob[0] - initial_prob[0]) + abs(new_initial_prob[1] - initial_prob[1]);
        long double muDiff = abs(new_mean[0] - mean[0]) + abs(new_mean[1] - mean[1]);
        long double sigmaDiff = abs(new_variance[0] - variance[0]) + abs(new_variance[1] - variance[1]);
        long double ADiff = 0.0;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                ADiff += abs(new_trans_matrix[i][j] - trans_matrix[i][j]);
            }
        }

        if((piDiff < tol) && (muDiff < tol) && (sigmaDiff < tol) && (ADiff < tol)){
            converged = true;
        }

        variance = new_variance;
        mean = new_mean;
        initial_prob = new_initial_prob;
        trans_matrix = new_trans_matrix;

        iter++;
    }

    bool ho = (mean[0] > mean[1]);

    vector<vector<long double>> alpha = forward_prob(trans_matrix, initial_prob, T, observation, mean, variance);
    vector<vector<long double>> beta = backward_prob(trans_matrix, initial_prob, T, observation, mean, variance);
    vector<vector<long double>> state_posterior = state_posterior_prob(alpha, T, beta);

    for (const auto& state : state_posterior) {
        if (!ho) {
            if (state[0] < state[1]) {
                cout << "Bull" << endl;
            }
            else {
                cout << "Bear" << endl;
            }
        }
        else {
            if (state[0] > state[1]) {
                cout << "Bull" << endl;
            }
            else {
                cout << "Bear" << endl;
            }
        }
    }

    // auto end = std::chrono::high_resolution_clock::now();
    // // Calculate the duration
    // std::chrono::duration<double> duration = end - start;

    // Output the duration
    // std::cout << "Elapsed time: " << duration.count() << " seconds" << std::endl;
    return 0;
}