#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <string.h>  
#include <time.h>

#define FLOAT_SIZE        32
#define DOUBLE_SIZE       64
#define LONG_DOUBLE_SIZE  128

void number_to_machine_code(void *num_ptr, size_t size, FILE *fp) {
    unsigned char *bytes = (unsigned char *)num_ptr;
    for (int byte = size - 1; byte >= 0; byte--) {
        for (int bit = 7; bit >= 0; bit--) {
            fprintf(fp, "%d", (bytes[byte] >> bit) & 1);
        }
        fprintf(fp, " ");
    }
}

void machine_code_to_number(void *bits_ptr, size_t size, void *result_ptr) {
    memcpy(result_ptr, bits_ptr, size);
}

double random_double(double min, double max) {
    double scale = rand() / (double)RAND_MAX;
    return min + scale * (max - min);
}


long double random_long_double(long double min, long double max) {
    long double scale = (long double)rand() / (long double)RAND_MAX;
    // добавляем extra случайные биты в младшие разряды
    long double extra = (long double)rand() / ((long double)RAND_MAX * (long double)RAND_MAX);
    return min + (scale + extra) * (max - min);
}

void generate_variants(int N, int K, int bits_type, int precision, float range_a, float range_b) {
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
            long double precise_num = random_long_double((long double)range_a, (long double)range_b);

            if (bits_type == FLOAT_SIZE) {
                float num = (float)precise_num;

                // 
                char buf[128];
                snprintf(buf, sizeof(buf), "%.*f", precision, num);
                float displayed;
                sscanf(buf, "%f", &displayed);
                double error = fabs((double)num - (double)displayed);

                fprintf(f_task,  "| %d | %.*f |\n", j, precision, num);
                fprintf(f_check, "| %d | %.*f | `", j, precision, num);
                number_to_machine_code(&num, sizeof(float), f_check);
                fprintf(f_check, "` | %.3e |\n", error);

            } else if (bits_type == DOUBLE_SIZE) {
                double num = (double)precise_num;

                // 
                char buf[128];
                snprintf(buf, sizeof(buf), "%.*f", precision, num);
                double displayed;
                sscanf(buf, "%lf", &displayed);
                long double error = fabsl((long double)num - (long double)displayed);

                fprintf(f_task,  "| %d | %.*f |\n", j, precision, num);
                fprintf(f_check, "| %d | %.*f | `", j, precision, num);
                number_to_machine_code(&num, sizeof(double), f_check);
                fprintf(f_check, "` | %.3Le |\n", error);

            } else { // long double (80-bit x86, хранится в 16 байтах с паддингом)
                long double num = precise_num;
                long double restored;
                machine_code_to_number(&num, sizeof(long double), &restored);

                // 
                char buf[128];
                snprintf(buf, sizeof(buf), "%.*Lf", precision, num);
                long double displayed;
                sscanf(buf, "%Lf", &displayed);
                long double error = fabsl(num - displayed);

                fprintf(f_task,  "| %d | %.*Lf |\n", j, precision, num);
                fprintf(f_check, "| %d | %.*Lf | `", j, precision, num);
                number_to_machine_code(&num, sizeof(long double), f_check);
                fprintf(f_check, "` | %.3Le |\n", error);
            }
        }

        fclose(f_task);
        fclose(f_check);
    }
}

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