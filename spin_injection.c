#include <stdio.h>
#include <math.h>

#define NN_HM (80)
#define NN_FM (50)
#define NN_NM (30)
#define NN_FM2 (60)
#define NN_NM2 (80)

int tridag(int N, double in[], double out[]);

int main(void)
{

  int i, j;

  double sigma_HM = 1.0e7; // プラチナ

  double sigma_FM = 1.0e7; // 鉄

  double sigma_NM = 6.0e7; // 銅

  double sigma_FM2 = 1.0e7; // 鉄

  double sigma_NM2 = 6.0e7; // 銅

  double kBT = (1.38e-23) * 300.0;

  double q = -1.602176634e-19;

  double theta = 0.1;

  double Jc = 1.0e10;

  double n = 6.0e23;

  double D = sigma_HM * kBT / (pow(q, 2.0) * n);

  double lambda_HM = 1.0e-9; // プラチナ

  double lambda_FM = 2.0e-9; // 鉄

  double lambda_NM = 200.0e-9; // 銅

  double lambda_FM2 = 2.0e-9; // 鉄

  double lambda_NM2 = 200.0e-9; // 銅

  double dz = 0.1e-9;

  double dt = 0.001;

  double s_HM[NN_HM + 2];

  double old_s_HM[NN_HM + 2];

  double new_s_HM[NN_HM + 2];

  double s_FM[NN_FM + 2];

  double old_s_FM[NN_FM + 2];

  double new_s_FM[NN_FM + 2];

  double s_NM[NN_NM + 2];

  double old_s_NM[NN_NM + 2];

  double new_s_NM[NN_NM + 2];

  double s_FM2[NN_FM2 + 2];

  double old_s_FM2[NN_FM2 + 2];

  double new_s_FM2[NN_FM2 + 2];

  double s_NM2[NN_NM2 + 2];

  double old_s_NM2[NN_NM2 + 2];

  double new_s_NM2[NN_NM2 + 2];

  double js_HM[NN_HM + 2];

  double js_FM[NN_FM + 2];

  double js_NM[NN_NM + 2];

  double js_FM2[NN_FM2 + 2];

  double js_NM2[NN_NM2 + 2];

  // HM 初期化設定

  for (i = 0; i <= NN_HM + 1; i++)
  {

    old_s_HM[i] = 0;

    old_s_HM[NN_HM + 1] = new_s_FM[0];

    js_HM[i] = 0;

    js_HM[NN_HM] = 2 * js_HM[NN_HM - 1] - js_HM[NN_HM - 2];
  }

  // FM 初期化設定

  for (i = 0; i <= NN_FM + 1; i++)
  {

    old_s_FM[i] = 0;

    new_s_FM[0] = old_s_HM[NN_HM + 1];

    js_FM[i] = 0;

    js_FM[0] = 2 * js_FM[1] - js_FM[2];
  }

  // NM 初期化設定

  for (i = 0; i <= NN_NM + 1; i++)
  {

    old_s_NM[i] = 0;

    new_s_NM[0] = old_s_FM[NN_FM + 1];

    js_NM[i] = 0;

    js_NM[0] = 2 * js_NM[1] - js_NM[2];
  }

  // FM2 初期化設定

  for (i = 0; i <= NN_FM2 + 1; i++)
  {

    old_s_FM2[i] = 0;

    new_s_FM2[0] = old_s_NM[NN_NM + 1];

    js_FM2[i] = 0;

    js_FM2[0] = 2 * js_FM2[1] - js_FM2[2];
  }

  // NM2 初期化設定

  for (i = 0; i <= NN_NM2 + 1; i++)
  {

    old_s_NM2[i] = 0;

    new_s_NM2[0] = old_s_NM2[NN_NM2 + 1];

    js_NM2[i] = 0;

    js_NM2[0] = 2 * js_NM2[1] - js_NM2[2];
  }

  // メインループ

  for (j = 1; j <= 10000000; j++)
  {

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ HM part start ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    // tridagの為の左端と右端の値代入（スピン流境界条件）

    s_HM[0] = ((sigma_HM * kBT * old_s_HM[1]) + (q * dz * theta * Jc)) / (sigma_HM * kBT);

    s_HM[NN_HM + 1] = old_s_HM[NN_HM] - ((q * dz) / (sigma_HM * kBT)) * (js_FM[0] + theta * Jc);

    // tridagの為の左端と右端の値代入（スピン蓄積密度境界条件）

    if (j / 1000 * 1000 == j)
    {

      s_HM[NN_HM + 1] = old_s_FM[1];
    }

    // tridagの為の値代入

    for (i = 1; i <= NN_HM; i++)
    {

      s_HM[i] = old_s_HM[i] * dz * dz / (lambda_HM * lambda_HM);
    }

    // tridagを用いてs（new_s_HM）を計算

    s_HM[1] = s_HM[1] - s_HM[0];

    s_HM[NN_HM] = s_HM[NN_HM] - s_HM[NN_HM + 1];

    s_HM[0] = 0;
    s_HM[NN_HM + 1] = 0;

    tridag(NN_HM, s_HM, new_s_HM);

    // new_s_HMを用いてold_s_HMを計算

    for (i = 1; i <= NN_HM; i++)
    {

      old_s_HM[i] = old_s_HM[i] + (new_s_HM[i] - old_s_HM[i]) * dt;
    }

    // 新たなold_s_HMからスピン流を計算

    for (i = 1; i <= NN_HM - 1; i++)
    {

      js_HM[i] = -((sigma_HM * kBT) / q) * ((old_s_HM[i + 1] - old_s_HM[i]) / dz) - theta * Jc;
    }

    // 補完

    js_HM[0] = 2 * js_HM[1] - js_HM[2];

    js_HM[NN_HM] = 2 * js_HM[NN_HM - 1] - js_HM[NN_HM - 2];

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ HM part finish ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ FM part start ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    s_FM[0] = ((q * dz) / (sigma_FM * kBT)) * js_HM[NN_HM] + old_s_FM[1];

    s_FM[NN_FM + 1] = old_s_FM[NN_FM] - ((q * dz) / (sigma_FM * kBT)) * (js_NM[0]);

    if (j / 100 * 100 == j)
    {

      s_FM[0] = old_s_HM[NN_HM];

      s_FM[NN_FM + 1] = old_s_NM[1];
    }

    for (i = 1; i <= NN_FM; i++)
    {

      s_FM[i] = old_s_FM[i] * dz * dz / (lambda_FM * lambda_FM);
    }

    s_FM[1] = s_FM[1] - s_FM[0];

    s_FM[NN_FM] = s_FM[NN_FM] - s_FM[NN_FM + 1];

    s_FM[0] = 0;
    s_FM[NN_FM + 1] = 0;

    tridag(NN_FM, s_FM, new_s_FM);

    for (i = 1; i <= NN_FM; i++)
    {

      old_s_FM[i] = old_s_FM[i] + (new_s_FM[i] - old_s_FM[i]) * dt;
    }

    for (i = 1; i <= NN_FM - 1; i++)
    {

      js_FM[i] = -((sigma_FM * kBT) / q) * ((old_s_FM[i + 1] - old_s_FM[i]) / dz);
    }

    js_FM[0] = 2 * js_FM[1] - js_FM[2];

    js_FM[NN_FM] = 2 * js_FM[NN_FM - 1] - js_FM[NN_FM - 2];

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ FM part finish ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ NM part start ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    s_NM[0] = ((q * dz) / (sigma_NM * kBT)) * js_FM[NN_FM] + old_s_NM[1];

    s_NM[NN_NM + 1] = old_s_NM[NN_NM] - ((q * dz) / (sigma_NM * kBT)) * (js_FM2[0]);

    if (j / 100 * 100 == j)
    {

      s_NM[0] = old_s_FM[NN_FM];

      s_NM[NN_NM + 1] = old_s_FM2[1];
    }

    for (i = 1; i <= NN_NM; i++)
    {

      s_NM[i] = old_s_NM[i] * dz * dz / (lambda_NM * lambda_NM);
    }

    s_NM[1] = s_NM[1] - s_NM[0];

    s_NM[NN_NM] = s_NM[NN_NM] - s_NM[NN_NM + 1];

    s_NM[0] = 0;
    s_NM[NN_NM + 1] = 0;

    tridag(NN_NM, s_NM, new_s_NM);

    for (i = 1; i <= NN_NM; i++)
    {

      old_s_NM[i] = old_s_NM[i] + (new_s_NM[i] - old_s_NM[i]) * dt;
    }

    for (i = 1; i <= NN_NM - 1; i++)
    {

      js_NM[i] = -((sigma_NM * kBT) / q) * ((old_s_NM[i + 1] - old_s_NM[i]) / dz);
    }

    js_NM[0] = 2 * js_NM[1] - js_NM[2];

    js_NM[NN_NM] = 2 * js_NM[NN_NM - 1] - js_NM[NN_NM - 2];

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ NM part finish ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ FM2 part start ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    s_FM2[0] = ((q * dz) / (sigma_FM2 * kBT)) * js_NM[NN_NM] + old_s_FM2[1];

    s_FM2[NN_FM2 + 1] = old_s_FM2[NN_FM2] - ((q * dz) / (sigma_FM2 * kBT)) * (js_NM2[0]);

    if (j / 100 * 100 == j)
    {

      s_FM2[0] = old_s_NM[NN_NM];

      s_FM2[NN_FM2 + 1] = old_s_NM2[1];
    }

    for (i = 1; i <= NN_FM2; i++)
    {

      s_FM2[i] = old_s_FM2[i] * dz * dz / (lambda_FM2 * lambda_FM2);
    }

    s_FM2[1] = s_FM2[1] - s_FM2[0];

    s_FM2[NN_FM2] = s_FM2[NN_FM2] - s_FM2[NN_FM2 + 1];

    s_FM2[0] = 0;
    s_FM2[NN_FM2 + 1] = 0;

    tridag(NN_FM2, s_FM2, new_s_FM2);

    for (i = 1; i <= NN_FM2; i++)
    {

      old_s_FM2[i] = old_s_FM2[i] + (new_s_FM2[i] - old_s_FM2[i]) * dt;
    }

    for (i = 1; i <= NN_FM2 - 1; i++)
    {

      js_FM2[i] = -((sigma_FM2 * kBT) / q) * ((old_s_FM2[i + 1] - old_s_FM2[i]) / dz);
    }

    js_FM2[0] = 2 * js_FM2[1] - js_FM2[2];

    js_FM2[NN_FM2] = 2 * js_FM2[NN_FM2 - 1] - js_FM2[NN_FM2 - 2];

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ FM2 part finish ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ NM2 part start ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★

    s_NM2[0] = ((q * dz) / (sigma_NM2 * kBT)) * js_FM2[NN_FM2] + old_s_NM2[1];

    s_NM2[NN_NM2 + 1] = old_s_NM2[NN_NM2];

    if (j / 100 * 100 == j)
    {

      s_NM2[0] = old_s_FM2[NN_FM2];
    }

    for (i = 1; i <= NN_NM2; i++)
    {

      s_NM2[i] = old_s_NM2[i] * dz * dz / (lambda_NM2 * lambda_NM2);
    }

    s_NM2[1] = s_NM2[1] - s_NM2[0];

    s_NM2[NN_NM2] = s_NM2[NN_NM2] - s_NM2[NN_NM2 + 1];

    s_NM2[0] = 0;
    s_NM2[NN_NM2 + 1] = 0;

    tridag(NN_NM2, s_NM2, new_s_NM2);

    for (i = 1; i <= NN_NM2; i++)
    {

      old_s_NM2[i] = old_s_NM2[i] + (new_s_NM2[i] - old_s_NM2[i]) * dt;
    }

    for (i = 1; i <= NN_NM2 - 1; i++)
    {

      js_NM2[i] = -((sigma_NM2 * kBT) / q) * ((old_s_NM2[i + 1] - old_s_NM2[i]) / dz);
    }

    js_NM2[0] = 2 * js_NM2[1] - js_NM2[2];

    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    // ★★★★★★★★★★ NM2 part finish ★★★★★★★★★★
    // ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
  }

  // HM出力

  for (i = 1; i <= NN_HM; i++)
  {

    printf("%d %s %e %s %e \n", i, "s", old_s_HM[i], "Js", -js_HM[i]);
  }

  // HM出力

  // FM出力

  for (i = 1; i <= NN_FM; i++)
  {

    printf("%d %s %e %s %e \n", i + NN_HM, "s", old_s_FM[i], "Js", -js_FM[i]);
  }

  // FM出力

  // NM出力

  for (i = 1; i <= NN_NM; i++)
  {

    printf("%d %s %e %s %e \n", i + NN_HM + NN_FM, "s", old_s_NM[i], "Js", -js_NM[i]);
  }

  // NM出力

  // FM2出力

  for (i = 1; i <= NN_FM2; i++)
  {

    printf("%d %s %e %s %e \n", i + NN_HM + NN_FM + NN_NM, "s", old_s_FM2[i], "Js", -js_FM2[i]);
  }

  // FM2出力

  // NM2出力

  for (i = 1; i <= NN_NM2; i++)
  {

    printf("%d %s %e %s %e \n", i + NN_HM + NN_FM + NN_NM + NN_FM2, "s", old_s_NM2[i], "Js", -js_NM2[i]);
  }

  // NM2出力
}

int tridag(int N, double in[], double out[])

{

  double x[10000];

  double g[10000];

  int j;

  double a, b, c, bet;

  a = 1.0;

  b = -2.0;

  c = 1.0;

  bet = b;

  x[1] = in[1] / bet;

  for (j = 2; j <= N; j++)

  {

    g[j] = c / (bet);

    bet = b - a * g[j];

    x[j] = (in[j] - a * x[j - 1]) / bet;
  }

  for (j = N - 1; j >= 1; j--)

  {

    x[j] = x[j] - g[j + 1] * x[j + 1];
  }

  for (j = 1; j <= N; j++)

  {

    out[j] = x[j];
  }
}