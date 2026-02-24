#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <string.h>  
#include <time.h>

#define FLOAT_SIZE   32
#define DOUBLE_SIZE  64
#define QUAD_SIZE    128

// float
void float_to_binary(float num, int binary[]) {
    uint32_t bits = 0;
    memcpy(&bits, &num, sizeof(float));
    for (int i = 0; i < FLOAT_SIZE; i++)
        binary[i] = (bits >> i) & 1;
}

float binary_to_float(int binary[]) {
    uint32_t bits = 0;
    for (int i = 0; i < FLOAT_SIZE; i++)
        bits |= ((uint32_t)binary[i] << i);
    float result;
    memcpy(&result, &bits, sizeof(float));
    return result;
}

// double
void double_to_binary(double num, int binary[]) {
    uint64_t bits = 0;
    memcpy(&bits, &num, sizeof(double));
    for (int i = 0; i < DOUBLE_SIZE; i++)
        binary[i] = (bits >> i) & 1;
}

double binary_to_double(int binary[]) {
    uint64_t bits = 0;
    for (int i = 0; i < DOUBLE_SIZE; i++)
        bits |= ((uint64_t)binary[i] << i);
    double result;
    memcpy(&result, &bits, sizeof(double));
    return result;
}

// 128
// Конвертируем double в 128-bit IEEE 754
void double_to_quad_bytes(double num, unsigned char bytes[16]) {
    memset(bytes, 0, 16);

    if (num == 0.0) return;

    // Знак
    int sign = (num < 0) ? 1 : 0;
    if (sign) num = -num;

    // Экспонента и мантисса из double
    int exp2;
    double mantissa = frexp(num, &exp2); // num = mantissa * 2^exp2, 0.5 <= mantissa < 1
    // Приводим к виду 1.xxx * 2^e
    exp2 -= 1;           // теперь mantissa in [1, 2)
    mantissa *= 2.0;     // mantissa = 1.xxx

    // Bias для quad = 16383
    int biased_exp = exp2 + 16383;

    // Записываем биты мантиссы (112 бит, неявная 1 не хранится)
    mantissa -= 1.0; // убираем неявную единицу
    // Записываем 112 бит мантиссы в bytes[0..13] (младшие байты)
    for (int i = 111; i >= 0; i--) {
        mantissa *= 2.0;
        int bit = (int)mantissa;
        mantissa -= bit;
        int byte_idx = i / 8;
        int bit_idx  = i % 8;
        bytes[byte_idx] |= (bit << bit_idx);
    }

    // Экспонента: биты 112..126 (15 бит)
    for (int i = 0; i < 15; i++) {
        int bit = (biased_exp >> i) & 1;
        int abs_bit = 112 + i;
        bytes[abs_bit / 8] |= (bit << (abs_bit % 8));
    }

    // Знак: бит 127
    if (sign)
        bytes[15] |= 0x80;
}

void quad_bytes_to_binary(unsigned char bytes[16], int binary[]) {
    for (int i = 0; i < QUAD_SIZE; i++)
        binary[i] = (bytes[i / 8] >> (i % 8)) & 1;
}

// Восстанавливаем double из 128-bit quad байт
double quad_bytes_to_double(unsigned char bytes[16]) {
    int sign = (bytes[15] >> 7) & 1;

    int biased_exp = 0;
    for (int i = 0; i < 15; i++) {
        int abs_bit = 112 + i;
        biased_exp |= (((bytes[abs_bit / 8] >> (abs_bit % 8)) & 1) << i);
    }

    if (biased_exp == 0) return 0.0;

    int exp2 = biased_exp - 16383;

    double mantissa = 1.0;
    double bit_val = 0.5;
    for (int i = 111; i >= 0; i--) {
        int bit = (bytes[i / 8] >> (i % 8)) & 1;
        mantissa += bit * bit_val;
        bit_val *= 0.5;
    }

    double result = ldexp(mantissa, exp2);
    return sign ? -result : result;
}

// вывод бинарного представления
void fprint_binary(FILE *fp, int binary[], int size) {
    fprintf(fp, "%d ", binary[size - 1]);

    int exp_start, exp_end;
    if (size == FLOAT_SIZE) {
        exp_start = 30; exp_end = 23;
    } else if (size == DOUBLE_SIZE) {
        exp_start = 62; exp_end = 52;
    } else {
        exp_start = 126; exp_end = 112;
    }

    for (int i = exp_start; i >= exp_end; i--)
        fprintf(fp, "%d", binary[i]);
    fprintf(fp, " ");

    for (int i = exp_end - 1; i >= 0; i--)
        fprintf(fp, "%d", binary[i]);
}

//  случайные числа 
double random_double(double min, double max) {
    double scale = rand() / (double)RAND_MAX;
    return min + scale * (max - min);
}

//  генерация вариантов 
void generate_variants(int N, int K, int bits_type, int precision, float range_a, float range_b) {
    int float_bits[FLOAT_SIZE];
    int double_bits[DOUBLE_SIZE];
    int quad_bits[QUAD_SIZE];

    for (int i = 1; i <= N; i++) {
        char task_filename[50], check_filename[50];
        sprintf(task_filename, "../tasks/variant_%d_task.md", i);
        sprintf(check_filename, "../checks/variant_%d_check.md", i);

        FILE *f_task  = fopen(task_filename, "w");
        FILE *f_check = fopen(check_filename, "w");

        if (!f_task || !f_check) {
            printf("Error creating files for variant %d\n", i);
            continue;
        }

        fprintf(f_task, "# Вариант %d\n\n", i);
        fprintf(f_task, "| N | Вещ число |\n|---|---|\n");

        fprintf(f_check, "# Проверка варианта %d\n\n", i);
        fprintf(f_check, "| N | Вещ число | Машинное представление | Ошибка |\n");
        fprintf(f_check, "|---|---|---|---|\n");

        for (int j = 1; j <= K; j++) {
            double precise_num = random_double((double)range_a, (double)range_b);

            if (bits_type == FLOAT_SIZE) {
                float num_f = (float)precise_num;
                float_to_binary(num_f, float_bits);
                double error = fabs(precise_num - (double)num_f);

                fprintf(f_task,  "| %d | %.*lf |\n", j, precision, precise_num);
                fprintf(f_check, "| %d | %.*lf | `", j, precision, precise_num);
                fprint_binary(f_check, float_bits, FLOAT_SIZE);
                fprintf(f_check, "` | %.3e |\n", error);

            } else if (bits_type == DOUBLE_SIZE) {
                double_to_binary(precise_num, double_bits);
                double restored = binary_to_double(double_bits);
                double error = fabs(precise_num - restored);

                fprintf(f_task,  "| %d | %.*lf |\n", j, precision, precise_num);
                fprintf(f_check, "| %d | %.*lf | `", j, precision, precise_num);
                fprint_binary(f_check, double_bits, DOUBLE_SIZE);
                fprintf(f_check, "` | %.3e |\n", error);

            } else { // 128 bit
                unsigned char quad_bytes[16];
                double_to_quad_bytes(precise_num, quad_bytes);
                quad_bytes_to_binary(quad_bytes, quad_bits);
                double restored = quad_bytes_to_double(quad_bytes);
                double error = fabs(precise_num - restored);

                fprintf(f_task,  "| %d | %.*lf |\n", j, precision, precise_num);
                fprintf(f_check, "| %d | %.*lf | `", j, precision, precise_num);
                fprint_binary(f_check, quad_bits, QUAD_SIZE);
                fprintf(f_check, "` | %.3e |\n", error);
            }
        }

        fclose(f_task);
        fclose(f_check);
    }
}

//  ввод данных 
void input_data(int *N, int *K, int *bits_type, int *precision, float *range_a, float *range_b) {
    FILE *inputFile = fopen("../input.txt", "r");
    if (!inputFile) {
        printf("input.txt не найден\n");
        return;
    }

    fscanf(inputFile, "%d", N);
    fscanf(inputFile, "%d", K);
    fscanf(inputFile, "%d", bits_type);
    fscanf(inputFile, "%f %f", range_a, range_b);
    fscanf(inputFile, "%d", precision);
    fclose(inputFile);

    printf("Кол-во вариантов: %d\n", *N);
    printf("Кол-во заданий: %d\n", *K);
    printf("Разрядность: %d\n", *bits_type);
    printf("Диапазон: [%f ; %f]\n", *range_a, *range_b);
    printf("Знаков после запятой: %d\n", *precision);
}

//  main 
int main() {
    srand(time(NULL));
    int N, K, bits_type, precision;
    float range_a, range_b;
    mkdir("../tasks",  0777);
    mkdir("../checks", 0777);
    input_data(&N, &K, &bits_type, &precision, &range_a, &range_b);
    generate_variants(N, K, bits_type, precision, range_a, range_b);
    printf("Готово!\n");
    return 0;
}