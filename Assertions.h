#define verify(expression)                                         \
    if (!expression){                                                \
        printf("\n%s%s", "Verification Failed: ", #expression);       \
        Exit();                                                        \
    }                                                                   \
