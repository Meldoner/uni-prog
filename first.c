#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <string.h>  
#include <time.h>

#define FLOAT_SIZE 32
#define DOUBLE_SIZE 64

void float_to_binary(float num, int binary[]) {
    for (int i = 0; i < FLOAT_SIZE; i++) {
        binary[i] = 0;
    }
    unsigned int* ptr = (unsigned int*)&num;
    unsigned int bits = *ptr;
    for (int i = 0; i < FLOAT_SIZE; i++) {
        binary[i] = (bits >> i) & 1;
    }
}

void double_to_binary(double num, int binary[]) {
    for (int i = 0; i < DOUBLE_SIZE; i++) {
        binary[i] = 0;
    }
    unsigned long long* ptr = (unsigned long long*)&num;
    unsigned long long bits = *ptr;
    for (int i = 0; i < DOUBLE_SIZE; i++) {
        binary[i] = (bits >> i) & 1;
    }
}

float binary_to_float(int binary[]) {
    unsigned int bits = 0;
    for (int i = 0; i < FLOAT_SIZE; i++) {
        bits |= (binary[i] << i);
    }
    float* ptr = (float*)&bits;
    return *ptr;
}

double binary_to_double(int binary[]) {
    unsigned long long bits = 0;
    for (int i = 0; i < DOUBLE_SIZE; i++) {
        bits |= ((unsigned long long)binary[i] << i);
    }
    double* ptr = (double*)&bits;
    return *ptr;
}

void fprint_binary(FILE *fp, int binary[], int size) {
    // 1. Знак (последний элемент массива)
    fprintf(fp, "%d ", binary[size - 1]);
    
    // Определяем границы для экспоненты
    int exp_start = (size == 32) ? 30 : 62;
    int exp_end   = (size == 32) ? 23 : 52;
    
    // 2. Экспонента
    for (int i = exp_start; i >= exp_end; i--) {
        fprintf(fp, "%d", binary[i]);
    }
    fprintf(fp, " "); // Пробел после экспоненты
    
    // 3. Мантисса
    for (int i = exp_end - 1; i >= 0; i--) {
        fprintf(fp, "%d", binary[i]);
    }
}

float random_float(float min, float max) {
    float scale = rand() / (float) RAND_MAX;
    return min + scale * ( max - min );
}

double random_double(double min, double max) {
    double scale = rand() / (double) RAND_MAX;
    return min + scale * ( max - min );
}

void generate_variants(int N, int K, int bits_type, int precision, float range_a, float range_b) {
    int float_bits[FLOAT_SIZE];
    int double_bits[DOUBLE_SIZE];
    for (int i = 1; i <= N; i++) {
        char task_filename[50], check_filename[50];
        sprintf(task_filename, "../tasks/variant_%d_task.md", i);
        sprintf(check_filename, "../checks/variant_%d_check.md", i);

        FILE *f_task = fopen(task_filename, "w");
        FILE *f_check = fopen(check_filename, "w");

        if (!f_task || !f_check) {
            printf("Error creating files for variant %d\n", i);
            continue;
        }

        // заголовки Markdown
        fprintf(f_task, "# Вариант %d\n\n", i);
        fprintf(f_task, "| N | Вещ число |\n");
        fprintf(f_task, "|---|---|\n");

        fprintf(f_check, "# Проверка варианта %d\n\n", i);
        fprintf(f_check, "| N | Вещ число | Машинное представление | Ошибка |\n");
        fprintf(f_check, "|---|---|---|---|\n");

        // 4. Генерация задач внутри варианта
        for (int j = 1; j <= K; j++) {
            double precise_num = random_double((double)range_a, (double)range_b);
            
            if (bits_type == 32) {
                float num_f = (float)precise_num; // Округление до 32 бит
                float_to_binary(num_f, float_bits);
                double error = fabs(precise_num - (double)num_f);

                // запись в файл задания
                fprintf(f_task, "| %d | %.*lf |\n", j, precision, precise_num);
                // запись в файл проверки
                fprintf(f_check, "| %d | %.*lf | `", j, precision, precise_num);
                fprint_binary(f_check, float_bits, FLOAT_SIZE);
                fprintf(f_check, "` | %.3e |\n", error);

            } else { // 64 bits
                double_to_binary(precise_num, double_bits);
                double restored = binary_to_double(double_bits);
                double error = precise_num - restored;
                
                // запись в файл задания
                fprintf(f_task, "| %d | %.*lf |\n", j, precision, precise_num);
                // запись в файл проверки
                fprintf(f_check, "| %d | %.*lf | `", j, precision, precise_num);
                fprint_binary(f_check, double_bits, DOUBLE_SIZE);
                fprintf(f_check, "` | %.3e |\n", error); // у double всегда ошибка = 0
            }
        }

        fclose(f_task);
        fclose(f_check);
    }
}

void input_data(int *N, int *K, int *bits_type, int *precision, float *range_a, float *range_b) {
    FILE *inputFile; 
    if ((inputFile = fopen("../input.txt", "r")) == NULL) {
        printf("input.txt не найден");
        return;
    }

    fscanf(inputFile, "%d", N);
    fscanf(inputFile, "%d", K);
    fscanf(inputFile, "%d", bits_type);
    fscanf(inputFile, "%f %f", range_a, range_b);
    fscanf(inputFile, "%d", precision);
    fclose(inputFile);

    printf("Кол-во вариантов: %d\n", *N);
    printf("Кол-во заданий на одного студента: %d\n", *K);
    printf("Разрядность: %d\n", *bits_type);
    printf("Диапазон: [%f ; %f]\n", *range_a, *range_b);
    printf("Кол-во знаков после запятой: %d\n", *precision);
}

int main() {
    srand(time(NULL));
    int N, K, bits_type, precision;
    float range_a, range_b;
    mkdir("tasks", 0777);
    mkdir("checks", 0777);
    input_data(&N, &K, &bits_type, &precision, &range_a, &range_b);
    generate_variants(N, K, bits_type, precision, range_a, range_b);

    printf("Готово!\n");
    return 0;
}