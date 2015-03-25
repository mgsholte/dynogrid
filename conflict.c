#include <stdio.h>

// a silly main function
int main(int argc, char *argv[]) {
	int i,b; // rename one of these

	// set a to be the average of the 1st 10 numbers
	b = 0;
	for(i = 0; i < 10; ++i) {
		b += i;
	}

	printf("the average of the 1st 10 numbers is: %d", b/10);
}
