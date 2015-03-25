#include <stdio.h>

// a silly main function
int main(int argc, char *argv[]) {
	int i,a; // rename one of these

	// set a to be the average of the 1st 10 numbers
	a = 0;
	for(i = 0; i < 10; ++i) {
		a += i;
	}

	printf("the average of the 1st 10 numbers is: %d", a/10);
}
