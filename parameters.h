#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <cmath>
#include <vector>
#include <C:\Users\Akira Ito\eigen-3.4.0\Eigen\Dense>

// シミュレーションのパラメータをまとめた構造体
struct SimulationParameters
{
    int K_ = 1405;                          // サブキャリア数 MODE1
    int L_ = 204;                          // 1フレームのシンボル数 L
    const int Q_ = 16;                           // 伝送路のインパルス応答のパス数 Q
    const double T_ = 252 * std::pow(10, -6);    // 有効シンボル長 T
    const double Tgi_ = 31.5 * std::pow(10, -6); // ガードインターバル長(1/8) Tgi
    const double Ts_ = T_ + Tgi_;               // シンボル全体の長さ
    const int NUMBER_OF_FFT = 2048;               // FFTポイント数
    const int NUMBER_OF_PILOT = 1;              // パイロットシンボル個数
    int NUMBER_OF_BIT;                          // 変調方式のビット数 (入力で設定)
    int NUMBER_OF_SYMBOLS;                      // シンボル数 (2^NUMBER_OF_BIT)
    int seed = 100;

    // パスの有無を制御するマスク
    std::vector<int> pathMask = {1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // DFT行列を生成する共通ロジック
    static Eigen::MatrixXcd generateW(int K, int Q, int FFT_size) {
        Eigen::MatrixXcd W(K, Q);
        for (int q = 0; q < Q; ++q) {
            for (int k = 0; k < K / 2; ++k) {
                // -26から-1番目のキャリヤ 
                W(k, q) = std::polar(1.0, -2.0 * M_PI * ((double)k - (double)K / 2.0) * (double)q / (double)FFT_size);
                // 1から26番目のキャリヤ 
                W(k + K / 2, q) = std::polar(1.0, -2.0 * M_PI * ((double)k + 1.0) * (double)q / (double)FFT_size);
            }
        }
        return W;
    }

    // シンボル番号 l, キャリア番号 k を入力とし、それが SP かどうかを判定する
    bool isScatteredPilot(int l, int k) const {
        // ISDB-T モード1 の SP 配置規則
        // 12キャリアおき、かつ4シンボル周期で3キャリアずつずれる
        return (k % 12 == 3 * (l % 4));
    }
};

// (k % 12 == 3 * (l % 4))

#endif /* PARAMETERS_H */