#include <stdio.h>
#include <flint/flint.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>

int main()
{
    // 設定參數 (你可調整此處)
    const slong n = 4; // 必須是2的次方
    slong poly_coeffs[] = {0, 2, 3, 1}; //  (由低次到高次)
    slong len = sizeof(poly_coeffs) / sizeof(poly_coeffs[0]);

    if (len > n) {
        printf("Error: 多項式次數不得超過 n-1。\n");
        return 1;
    }

    // 初始化數域 context (Q[x]/(x^n+1))
    fmpq_poly_t minimal_poly;
    fmpq_poly_init(minimal_poly);
    fmpq_poly_set_coeff_si(minimal_poly, n, 1); // x^n
    fmpq_poly_set_coeff_si(minimal_poly, 0, 1); // +1 (x^n + 1)

    nf_t nf;
    nf_init(nf, minimal_poly);

    // 設定數域元素（待計算的多項式）
    nf_elem_t elem;
    nf_elem_init(elem, nf);

    fmpq_poly_t poly_input;
    fmpq_poly_init(poly_input);

    for (slong i = 0; i < len; ++i) {
        fmpq_poly_set_coeff_si(poly_input, i, poly_coeffs[i]);
    }

    nf_elem_set_fmpq_poly(elem, poly_input, nf);

    // 計算範數
    fmpq_t norm;
    fmpq_init(norm);
    nf_elem_norm(norm, elem, nf);

    // 輸出結果
    printf("Norm of polynomial: ");
    fmpq_poly_print_pretty(poly_input,"xbar");
    printf("\n");
    
    printf("is: ");
    fmpq_print(norm);
    printf("\n");

    // 清理資源
    fmpq_clear(norm);
    nf_elem_clear(elem, nf);
    nf_clear(nf);
    fmpq_poly_clear(poly_input);
    fmpq_poly_clear(minimal_poly);

    return 0;
}

