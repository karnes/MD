/*
 *	This is the main driver program for the molecular dynamics program.
 *	The user can specify a script of actions for the program to undertake,
 *	and this main program reads the script and sends the actions to 
 *	appropriate subroutines.
 */

#include <md.h>

main(argc,argv)
	int	argc;
	char	*argv[];
{
	FILE	*fp;
	char	buf[132]; /* moved action[5] to md.h */
	if	(argc != 2) {
		fprintf(stderr,"usage: %s  <setup file>\n", argv[0]);
		exit(0);
	}
	
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "main: cannot open %s\n", argv[1]);
		exit(1);
	}
	
/*** Search through the setup file for the action word and call the routine***/

	while (fgets(buf, 132, fp) != NULL) {
		sscanf(buf, "%s", action);
		if (strcmp(action, ".MIN.") == 0) {
			doMin(fp); /* do minimization */
		}
		else if (strcmp(action, ".EQU.") == 0) {
			doEqu(fp); /* do equilibriation */
		}
		else if (strcmp(action, ".RUN.") == 0) {
			doTraj(fp); /* run trajectory using equil. inp. files */
		}
		else if (strcmp(action, ".STR.") == 0) {
			reStart(fp); /* restart trajectory */
		}
		else if (strcmp(action, ".SPL.") == 0) {
			doSpec(fp); /* special run */
		}
	}
	
	fclose(fp);
	
}
