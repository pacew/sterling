#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <proj_api.h>

#define MASSGIS_URL "http://wsgw.mass.gov/data/gispub/images/coq2008_30cm_sid"

FILE *runf;
FILE *getf;

char output_dir[1000];
char tiffs_dir[1000];

projPJ pj_condor;
projPJ pj_massgis;

void make_tile (int tile_x, int tile_y);

/* massgis is serving these bad files ... ignore them */
char *bad_files[] = {
	"17279380",
	"17429380",
	"16228735",
	"19679080",
	"19378930",
};

char *hires_tiles[] = {
};

int
is_hires_tile (int hpos, int vpos)
{
	char buf[100];
	int idx;

	sprintf (buf, "%02d%02d", hpos, vpos);
	for (idx = 0; idx < sizeof hires_tiles / sizeof hires_tiles[0]; idx++) {
		if (strcmp (hires_tiles[idx], buf) == 0)
			return (1);
	}
	return (0);
}

int
should_ignore_file (char *name)
{
	char buf[1000];
	char *p;
	int idx;

	strncpy (buf, name, sizeof buf);
	buf[sizeof buf-1] = 0;
	if ((p = strchr (buf, '.')) != NULL)
		*p = 0;

	for (idx = 0; idx < sizeof bad_files / sizeof bad_files[0]; idx++) {
		if (strcmp (buf, bad_files[idx]) == 0)
			return (1);
	}
	return (0);
}

void
usage (void)
{
	fprintf (stderr, "usage: mktiles [tilenum]...\n");
	fprintf (stderr, "  tilenum is 1208, etc\n");
	exit (1);
}

void
pave_path (char *filename)
{
	char prefix[1000];
	char *p;
	int c;

	if (strlen (filename) > 900) {
		fprintf (stderr, "filename too long\n");
		exit (1);
	}

	strcpy (prefix, filename);

	p = prefix;
	while (*p == '/')
		p++;
	while (1) {
		while (*p && *p != '/')
			p++;
		if (*p == 0)
			break;
		c = *p;
		*p = 0;
		mkdir (prefix, 0777);
		*p = c;
		while (*p == '/')
			p++;
	}
}

char *hdr_map_projection;
char *hdr_ellipsoid;
char *hdr_left_map_x;
char *hdr_lower_map_y;
char *hdr_right_map_x;
char *hdr_upper_map_y;
char *hdr_number_of_rows;
char *hdr_number_of_columns;

struct name_valp {
	char *name;
	char **valp;
};

struct name_valp hdr[] = {
	{ "map_projection", &hdr_map_projection },
	{ "ellipsoid", &hdr_ellipsoid },
	{ "left_map_x", &hdr_left_map_x },
	{ "lower_map_y", &hdr_lower_map_y },
	{ "right_map_x", &hdr_right_map_x },
	{ "upper_map_y", &hdr_upper_map_y },
	{ "number_of_rows", &hdr_number_of_rows },
	{ "number_of_columns", &hdr_number_of_columns },
	{ NULL, NULL },
};

void
save_val (struct name_valp *tbl, char *name, char *val)
{
	struct name_valp *np;

	for (np = tbl; np->name; np++) {
		if (strcmp (np->name, name) == 0) {
			if ((*np->valp = strdup (val)) == NULL) {
				fprintf (stderr, "out of memory\n");
				exit (1);
			}
		}
	}
}

int left_map_x;
int lower_map_y;
int right_map_x;
int upper_map_y;
int number_of_rows;
int number_of_columns;

void
read_3dem_hdr (char *filename)
{
	FILE *f;
	char buf[1000];
	int len;
	char *p;
	char *name, *val;
	int err;
	struct name_valp *np;

	if ((f = fopen (filename, "r")) == NULL) {
		fprintf (stderr, "can't open %s\n", filename);
		exit (1);
	}

	while (fgets (buf, sizeof buf, f) != NULL) {
		len = strlen (buf);
		while (len > 0 && isspace (buf[len-1]))
			buf[--len] = 0;

		p = buf;
		while (isspace (*p))
			p++;
		name = p;
		while (*p && ! isspace (*p))
			p++;
		if (*p)
			*p++ = 0;
		while (isspace (*p) || *p == '=')
			p++;
		val = p;
		
		save_val (hdr, name, val);
	}
		
	fclose (f);


	err = 0;
	for (np = hdr; np->name; np++) {
		if (np->valp == NULL) {
			printf ("%s: missing value for %s\n",
				filename, np->name);
			err = 1;
		}
	}
	if (err)
		exit (1);

	if (strcmp (hdr_map_projection, "UTM Zone 19N") != 0) {
		fprintf (stderr, "bad projection\n");
		exit (1);
	}
	if (strcmp (hdr_ellipsoid, "WGS84") != 0) {
		fprintf (stderr, "bad ellipsoid\n");
		exit (1);
	}

	left_map_x = atoi (hdr_left_map_x);
	lower_map_y = atoi (hdr_lower_map_y);
	right_map_x = atoi (hdr_right_map_x);
	upper_map_y = atoi (hdr_upper_map_y);

	number_of_rows = atoi (hdr_number_of_rows);
	number_of_columns = atoi (hdr_number_of_columns);
}

/* http://en.wikipedia.org/wiki/World_file */

struct worldfile {
	struct worldfile *next;
	char *name;
	char *raster_name;
	double left, top, right, bottom;
	int used;
};

struct worldfile *world_files;

char *massgis_dir;

void
read_world_files (void)
{
	DIR *dir;
	struct dirent *dp;
	char *p;
	char fullname[1000];
	FILE *inf;
	double a, b, c, d, e, f;
	struct worldfile *wp;
	char raster_name[1000];
	int pixels_wide, pixels_high;
	double pixel_width, pixel_height;

	if ((dir = opendir (massgis_dir)) == NULL) {
		fprintf (stderr, "can't open directory %s\n",
			massgis_dir);
		exit (1);
	}

	while ((dp = readdir (dir)) != NULL) {
		if ((p = strrchr (dp->d_name, '.')) == NULL)
			continue;
		if (strcmp (p, ".tfw") != 0)
			continue;

		snprintf (fullname, sizeof fullname,
			  "%s/%s", massgis_dir, dp->d_name);
		if ((inf = fopen (fullname, "r")) == NULL) {
			fprintf (stderr, "can't open %s\n", fullname);
			exit (1);
		}
		if (fscanf (inf, "%lf %lf %lf %lf %lf %lf",
			    &a, &d, &b, &e, &c, &f) != 6) {
			fprintf (stderr, "%s: bad format\n", fullname);
			exit (1);
		}

		if (d != 0 || b != 0
		    || (fabs (a) - fabs (e) > 1e-6)) {
			fprintf (stderr, "%s: bad world file\n", fullname);
			exit (1);
		}

		strcpy (raster_name, dp->d_name);
		if ((p = strrchr (raster_name, '.')) != NULL)
			*p = 0;
		strcat (raster_name, ".tif");

		wp = calloc (1, sizeof *wp);
		wp->raster_name = strdup (raster_name);
		pixel_width = a;
		pixel_height = e;
		wp->left = c;
		wp->top = f;
			 
		pixels_wide = 5000;
		pixels_high = 5000;

		wp->right = wp->left + pixels_wide * pixel_width;
		wp->bottom = wp->top + pixels_high * pixel_height;

		wp->next = world_files;
		world_files = wp;

		fclose (inf);
	}

	closedir (dir);
}

int tiles_wide;
int tiles_high;

#define HIRES_TILE_SIZE 4096
#define NORMAL_TILE_SIZE 2048
#define SMALL_TILE_SIZE 512

int tile_list[1000];
int tile_list_used;

void
read_imgcat (void)
{
	FILE *inf;
	char buf[1000];
	char *p;
	char filename[1000];
	double xmin, ymin, xmax, ymax;
	struct worldfile *wp;
	char imgcat[1000];

	/*
	 * IMGCAT_COQ2008_2009SID_30CM.dbf comes from massgis
	 * pgdbf from the pgdbf package converts it to this text form
	 */
	sprintf (imgcat, "%s/imgcat.txt", massgis_dir);
	if ((inf = fopen (imgcat, "r")) == NULL) {
		fprintf (stderr, "can't open %s\n", imgcat);
		exit (1);
	}
	while (fgets (buf, sizeof buf, inf) != NULL) {
		if (sscanf (buf, "%s %lf %lf %lf %lf",
			    filename, &xmin, &ymin, &xmax, &ymax) != 5)
			continue;
		
		if (should_ignore_file (filename))
			continue;

		wp = calloc (1, sizeof *wp);

		if ((p = strrchr (filename, '.')) != NULL)
			*p = 0;
		wp->name = strdup (filename);

		wp->raster_name = malloc (strlen (wp->name) + 100);
		sprintf (wp->raster_name, "%s.tif", wp->name);

		wp->left = xmin;
		wp->top = ymax;
		wp->right = xmax;
		wp->bottom = ymin;

		wp->next = world_files;
		world_files = wp;
		    
	}
}

struct get_list {
	struct get_list *next;
	struct worldfile *wp;
};

struct get_list *get_list;

void
add_to_get_list (struct worldfile *wp)
{
	struct get_list *gp;
	for (gp = get_list; gp; gp = gp->next) {
		if (strcmp (gp->wp->name, wp->name) == 0)
			return;
	}
	gp = calloc (1, sizeof *gp);
	gp->wp = wp;
	gp->next = get_list;
	get_list = gp;
}

int
main (int argc, char **argv)
{
	int c;
	int tile_hpos, tile_vpos;
	char *hdrfile;
	int idx;
	char curdir[1000];
	struct get_list *gp;

	while ((c = getopt (argc, argv, "")) != EOF) {
		switch (c) {
		default:
			usage ();
		}
	}

	tile_list_used = 0;

	if (optind < argc) {
		while (optind < argc) {
			tile_list[tile_list_used++]
				= strtol (argv[optind++], NULL, 10);
		}
	} else {
		int hpos, vpos;

		for (hpos = 6; hpos <= 11; hpos++) {
			for (vpos = 0; vpos <= 15; vpos++) {
				tile_list[tile_list_used++] = hpos * 100 + vpos;
			}
		}
	}

	if ((runf = fopen ("TMP.run", "w")) == NULL) {
		fprintf (stderr, "can't create TMP.run\n");
		exit (1);
	}

	if ((getf = fopen ("TMP.get", "w")) == NULL) {
		fprintf (stderr, "can't create TMP.get\n");
		exit (1);
	}

	getcwd (curdir, sizeof curdir);
	if (strlen (curdir) > 900) {
		fprintf (stderr, "crazy current directory\n");
		exit (1);
	}

	sprintf (tiffs_dir, "%s/tiffs", curdir);

	hdrfile = "Sterling.hdr";
	massgis_dir = "/home/pace/sterling-data/massgis";

	fprintf (getf, "cd %s\n", massgis_dir);

	read_3dem_hdr (hdrfile);

	read_imgcat ();

	if (number_of_rows % 64 != 0
	    || number_of_columns % 64 != 0) {
		printf ("%s: rows and cols must be multiples of 64\n",
			hdrfile);
		exit (1);
	}

	tiles_wide = number_of_columns / 64;
	tiles_high = number_of_rows / 64;

	pj_condor = pj_init_plus ("+proj=utm"
				  " +zone=19"
				  " +ellps=WGS84"
				  " +datum=WGS84"
				  " +units=m"
				  " +no_defs");

	/* epsg:26986 NAD83 / Massachusetts Mainland */
	pj_massgis = pj_init_plus ("+proj=lcc"
				   " +lat_1=42.68333333333333"
				   " +lat_2=41.71666666666667"
				   " +lat_0=41"
				   " +lon_0=-71.5"
				   " +x_0=200000"
				   " +y_0=750000"
				   " +ellps=GRS80"
				   " +datum=NAD83"
				   " +units=m"
				   " +no_defs");

	for (tile_hpos = 0; tile_hpos < tiles_wide; tile_hpos++) {
		for (tile_vpos = 0; tile_vpos < tiles_high; tile_vpos++) {
			if (tile_list_used) {
				int tile;

				tile = tile_hpos * 100 + tile_vpos;

				for (idx = 0; idx < tile_list_used; idx++) {
					if (tile_list[idx] == tile)
						break;
				}
				if (idx == tile_list_used)
					continue;
			}
			make_tile (tile_hpos, tile_vpos);
		}
	}

	for (gp = get_list; gp; gp = gp->next) {
		struct worldfile *wp;

		wp = gp->wp;
		fprintf (getf, "wget '%s/%s.zip'\n",
			 MASSGIS_URL, wp->name);
		fprintf (getf, "./do-mrsid %s.zip\n", wp->name);
	}

	return (0);
}

void
make_tile (int tile_hpos, int tile_vpos)
{
	double tile_width, tile_height;
	int tiles_from_left;
	double left_tile_x, right_tile_x;
	double upper_tile_y, lower_tile_y;
	projUV coords;
	double left_photo, top_photo, right_photo, bottom_photo;
	struct worldfile *wp;
	char *source_proj, *dest_proj;
	FILE *f;
	char *cmd;
	size_t cmdsize;
	char outname[1000];
	char tname[1000];
	int size;
	int need_get;
	char inname[1000];

	sprintf (tname, "t%02d%02d", tile_hpos, tile_vpos);

	tile_width = (double)(right_map_x - left_map_x) / tiles_wide;
	tiles_from_left = tiles_wide - tile_hpos - 1;
	left_tile_x = left_map_x + tiles_from_left * tile_width;
	right_tile_x = left_tile_x + tile_width;

	tile_height = (double)(upper_map_y - lower_map_y) / tiles_high;
	lower_tile_y = lower_map_y + tile_vpos * tile_height;
	upper_tile_y = lower_tile_y + tile_height;

	coords.u = left_tile_x;
	coords.v = lower_tile_y;
	coords = pj_fwd (pj_inv (coords, pj_condor), pj_massgis);
	left_photo = coords.u;
	bottom_photo = coords.v;

	coords.u = right_tile_x;
	coords.v = upper_tile_y;
	coords = pj_fwd (pj_inv (coords, pj_condor), pj_massgis);
	right_photo = coords.u;
	top_photo = coords.v;

	for (wp = world_files; wp; wp = wp->next)
		wp->used = 0;

	need_get = 0;
	for (wp = world_files; wp; wp = wp->next) {
		if (wp->right < left_photo)
			continue;
		if (wp->left > right_photo)
			continue;
		if (wp->top < bottom_photo)
			continue;
		if (wp->bottom > top_photo)
			continue;
		wp->used = 1;

		sprintf (inname, "%s/%s", massgis_dir, wp->raster_name);
		if (access (inname, R_OK) < 0) {
			need_get = 1;

			add_to_get_list (wp);
		}
	}

	if (need_get)
		return;

	source_proj = "EPSG:26986";
	dest_proj = "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84"
		" +units=m +no_defs";

	size = NORMAL_TILE_SIZE;
	if (is_hires_tile (tile_hpos, tile_vpos))
		size = HIRES_TILE_SIZE;

	f = open_memstream (&cmd, &cmdsize);
	fprintf (f, "(cd %s && ", massgis_dir);
	fprintf (f, " gdalwarp");
	fprintf (f, " -s_srs '%s'", source_proj);
	fprintf (f, " -t_srs '%s'", dest_proj);
	fprintf (f, " -ts %d %d", size, size);
	fprintf (f, " -te %.12g %.12g %.12g %.12g",
			 left_tile_x, lower_tile_y,
		 right_tile_x, upper_tile_y);
	for (wp = world_files; wp; wp = wp->next) {
		if (! wp->used)
			continue;
		fprintf (f, " %s", wp->raster_name);
	}

	sprintf (tname, "t%02d%02d", tile_hpos, tile_vpos);

	sprintf (outname, "%s/%s.tif", tiffs_dir, tname);
	pave_path (outname);

	fprintf (f, " %s", outname);
	fprintf (f, ")");
	fclose (f);

	fprintf (runf, "%s\n", cmd);

	fprintf (runf, "tifftopnm %s | pnmflip -rotate180 > TMP.ppm\n",
		 outname);

	sprintf (outname, "Textures/%s.bmp", tname);
	pave_path (outname);
	fprintf (runf, "ppmtobmp TMP.ppm > %s\n", outname);

	sprintf (outname, "Textures/Small/%s.bmp", tname);
	pave_path (outname);
	fprintf (runf, "pnmscale -xsize=%d TMP.ppm | ppmtobmp > %s\n",
		 SMALL_TILE_SIZE, outname);
	
	fprintf (runf, "\n");
}
