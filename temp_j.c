#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[]) {
	// // int i;  // loop index var
	// int x_divs, y_divs, z_divs;
	// char x_temp = argv[1];

	// x_divs = atoi(x_temp);
	// y_divs = atoi(argv[2]);
	// z_divs = atoi(argv[3]);


	// printf("y_divs is: %d\n", y_divs);

	// printf("z_divs is: %d\n", z_divs);

	int num;
	sscanf (argv[1],"%d",&num);
	printf("x_divs is: %d\n", num);


	return 0;
}//end main