#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

// OBJECT STRUCTURE THAT ALLOWS FOR ALL 3 OBJECTS
typedef struct {
  int kind; // 0 = camera, 1 = sphere, 2 = plane, 3 = light
  double ns;
  union {
    struct {
    	double width;
    	double height;
    } camera;
    struct {
		double diffuse_color[3];
		double specular_color[3];
     	double position[3];
     	double radius;
    } sphere;
    struct {
		double diffuse_color[3];
		double specular_color[3];
     	double position[3];
    	double normal[3];
    } plane;
	struct {
		double color[3];
		double theta;
		double radial_a2;
		double radial_a1;
		double radial_a0;
		double angular_a0;
		double position[3];
		double direction[3];
	} light;
  };
} Object;

// PIXEL STRUCTURE CREATED
typedef struct Pixel{
	unsigned char red;
	unsigned char green;
	unsigned char blue;
} Pixel;

// OBJECT ARRAY TO READ FROM JSON FILE INTO
Object* object_array[128];
Object* lights[128];
int obj = 0;
int light = 0;
int line = 1;

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);    
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }  
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);      
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);      
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  // Error check this..
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}


void read_scene(char* filename) {
  int c;
  FILE* json = fopen(filename, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }
  
  skip_ws(json);
  
  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects

  while (1) {
  	c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return;
    }
    if (c == '{') {
      skip_ws(json);
    
      // Parse the object
      char* key = next_string(json);
      if (strcmp(key, "type") != 0) {
		fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
		exit(1);
      }

      skip_ws(json);

      expect_c(json, ':');

      skip_ws(json);

      char* value = next_string(json);
  	  object_array[obj] = malloc(sizeof(Object));
	  lights[light] = malloc(sizeof(Object));
	  // IDENTIFYING OBJECT TYPES AND BUILDING OBJECT
      if (strcmp(value, "camera") == 0) {
				(*object_array[obj]).kind = 0;
				printf("Found camera\n");
      } else if (strcmp(value, "sphere") == 0) {
				(*object_array[obj]).kind = 1;
				printf("Found sphere\n");
      } else if (strcmp(value, "plane") == 0) {
				(*object_array[obj]).kind = 2;
				printf("Found plane\n");
      } else if (strcmp(value, "light") == 0) {
				(*object_array[obj]).kind = 3;
				(*lights[light]).kind = 3;
				printf("Found light\n");
      } else {
				fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
				exit(1);
      }

      skip_ws(json);

      while (1) {
		c = next_c(json);
		if (c == '}') {
	  		// stop parsing this object
	  		//object_array[obj] = &new;
			obj++;
			if((*object_array[obj-1]).kind == 3){
				light++;					
			}
	  		break;
	  	} else if (c == ',') {
	  		// read another field
			skip_ws(json);
			char* key = next_string(json);
		  	skip_ws(json);
		  	expect_c(json, ':');
		  	skip_ws(json);
			// BUILDING OBJECT DOUBLE FIELDS
	  		if ((strcmp(key, "width") == 0) ||
				(strcmp(key, "height") == 0) ||
				(strcmp(key, "radius") == 0) ||
				(strcmp(key, "theta") == 0) ||
				(strcmp(key, "radial-a0") == 0) ||
				(strcmp(key, "radial-a1") == 0) ||
				(strcmp(key, "radial-a2") == 0)) {
				double value = next_number(json);
				if(strcmp(key, "width") == 0){
					if((*object_array[obj]).kind == 0) (*object_array[obj]).camera.width = value;
				}
				else if(strcmp(key, "height") == 0){
					if((*object_array[obj]).kind == 0) (*object_array[obj]).camera.height = value;
				}
				else if(strcmp(key, "radius") == 0){
					(*object_array[obj]).sphere.radius = value;
				}
				else if(strcmp(key, "theta") == 0){
					if((*object_array[obj]).kind == 3) {
						(*object_array[obj]).light.theta = value;
						(*lights[light]).light.theta = value;
					}
				}
				else if(strcmp(key, "radial-a0") == 0){
					if((*object_array[obj]).kind == 3) {
						(*object_array[obj]).light.radial_a0 = value;
						(*lights[light]).light.radial_a0 = value;
					}
				}
				else if(strcmp(key, "radial-a1") == 0){
					if((*object_array[obj]).kind == 3) {
						(*object_array[obj]).light.radial_a1 = value;
						(*lights[light]).light.radial_a1 = value;
					}
				}
				else if(strcmp(key, "radial-a2") == 0){
					if((*object_array[obj]).kind == 3) {
						(*object_array[obj]).light.radial_a2 = value;
						(*lights[light]).light.radial_a2 = value;
					}
				}
				// BUILDING OBJECT VECTOR FIELDS
	  		} else if ((strcmp(key, "color") == 0) ||
		    	(strcmp(key, "position") == 0) ||
		    	(strcmp(key, "normal") == 0) ||
				(strcmp(key, "diffuse_color") == 0) ||
				(strcmp(key, "specular_color") == 0)) {
	    		double* value = next_vector(json);
				if(strcmp(key, "color") == 0){
					if((*object_array[obj]).kind == 3){
						(*object_array[obj]).light.color[0] = value[0];
						(*object_array[obj]).light.color[1] = value[1];
						(*object_array[obj]).light.color[2] = value[2];

						(*lights[light]).light.color[0] = value[0];
						(*lights[light]).light.color[1] = value[1];
						(*lights[light]).light.color[2] = value[2];
					}
				}
				else if(strcmp(key, "position") == 0){
					if((*object_array[obj]).kind == 1){
						(*object_array[obj]).sphere.position[0] = value[0];
						(*object_array[obj]).sphere.position[1] = value[1];
						(*object_array[obj]).sphere.position[2] = value[2];
					}
					else if((*object_array[obj]).kind == 2){
						(*object_array[obj]).plane.position[0] = value[0];
						(*object_array[obj]).plane.position[1] = value[1];
						(*object_array[obj]).plane.position[2] = value[2];
					}
					else if((*object_array[obj]).kind == 3){
						(*object_array[obj]).light.position[0] = value[0];
						(*object_array[obj]).light.position[1] = value[1];
						(*object_array[obj]).light.position[2] = value[2];

						(*lights[light]).light.position[0] = value[0];
						(*lights[light]).light.position[1] = value[1];
						(*lights[light]).light.position[2] = value[2];
					}
				}
				else if(strcmp(key, "normal") == 0){
					(*object_array[obj]).plane.normal[0] = value[0];
					(*object_array[obj]).plane.normal[1] = value[1];
					(*object_array[obj]).plane.normal[2] = value[2];
				}
				else if(strcmp(key, "diffuse_color") == 0){
					if((*object_array[obj]).kind == 1){
						(*object_array[obj]).sphere.diffuse_color[0] = value[0];
						(*object_array[obj]).sphere.diffuse_color[1] = value[1];
						(*object_array[obj]).sphere.diffuse_color[2] = value[2];
					}
					else if((*object_array[obj]).kind == 2){
						(*object_array[obj]).plane.diffuse_color[0] = value[0];
						(*object_array[obj]).plane.diffuse_color[1] = value[1];
						(*object_array[obj]).plane.diffuse_color[2] = value[2];
					}
				}
				else if(strcmp(key, "specular_color") == 0){
					if((*object_array[obj]).kind == 1){
						(*object_array[obj]).sphere.specular_color[0] = value[0];
						(*object_array[obj]).sphere.specular_color[1] = value[1];
						(*object_array[obj]).sphere.specular_color[2] = value[2];
					}
					else if((*object_array[obj]).kind == 2){
						(*object_array[obj]).plane.specular_color[0] = value[0];
						(*object_array[obj]).plane.specular_color[1] = value[1];
						(*object_array[obj]).plane.specular_color[2] = value[2];
					}
				}
				else if(strcmp(key, "direction") == 0){
					if((*object_array[obj]).kind == 3){
						(*object_array[obj]).light.direction[0] = value[0];
						(*object_array[obj]).light.direction[1] = value[1];
						(*object_array[obj]).light.direction[2] = value[2];

						(*lights[light]).light.direction[0] = value[0];
						(*lights[light]).light.direction[1] = value[1];
						(*lights[light]).light.direction[2] = value[2];
					}
				}
				// ERROR CHECK
	 		 } else {
	   			 fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
		   		 key, line);
	   			 //char* value = next_string(json);
	 		 }
	 		 (*object_array[obj]).ns = 20.0;
	 		 if ((*object_array[obj]).kind == 3 && (*object_array[obj]).light.theta == 0){
				(*object_array[obj]).light.theta = 0;
				(*lights[light]).light.theta = 0;
	  		}
	  		if ((*object_array[obj]).kind == 3 && (*object_array[obj]).light.angular_a0 == 0){
				(*object_array[obj]).light.angular_a0 = 0;
				(*lights[light]).light.angular_a0 = 0;
	 		 }
	 		 skip_ws(json);
		} else {
	 		 fprintf(stderr, "Error: Unexpected value on line %d\n", line);
	 		 exit(1);
		}
	}
    skip_ws(json);
    c = next_c(json);
    if (c == ',') {
		// noop
		skip_ws(json);
    } else if (c == ']') {
		fclose(json);
		object_array[obj] = NULL;
		return;
    } else {
		fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
		exit(1);
    }
    }
  }
}


///////////////////////////////////////////////////////////////
// BEGINNING OF RAYCASTING FUNCTION
///////////////////////////////////////////////////////////////

double sqr(double v) {
  return v*v;
}

double dot(double* a, double* b) {
	double result;
	result = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	return result;
}

double* diffuse(double* Kd, double* Il, double* N, double* L){
	double* diff;
	double dot1 = N[0]*L[0] + N[1]*L[1] + N[2]*L[2];
	if (dot > 0){
		diff[0] = dot1*Kd[0]*Il[0];
		diff[1] = dot1*Kd[1]*Il[1];
		diff[2] = dot1*Kd[2]*Il[2];
		return diff;
	}
	else {
		diff[0] = 0;
		diff[1] = 0;
		diff[2] = 0;
		return diff;
	}
}

double* specular(double* Ks, double* Il, double* V, double* R, double* N, double* L){
	double dot1 = dot(V,R);
	double dot2 = dot(N,L);
	double* result;
	if (dot1 > 0 && dot2 > 0) {
		result[0] = pow(dot1,20)*Ks[0]*Il[0];
		result[1] = pow(dot1,20)*Ks[1]*Il[1];
		result[2] = pow(dot1,20)*Ks[2]*Il[2];
		return result;
	}
	else {
		result[0] = 0;
		result[1] = 0;
		result[2] = 0;
		return result;
	}
}

double frad(double* r, double d){
	double determinant = ((r[0] * sqr(d)) + (r[1] * d) + r[2]);
	if (determinant == 0) {
		return 0;
	}
	else{
		return 1/determinant;
	}
}

double fang(double a, double t, double* Rdn, double* direction){
	if (direction[0] == 0 && direction[1] == 0 && direction[2] == 0){
		return 1;
	}
	double value = dot(Rdn, direction);
	if (value < cos(t)) {
		return 0;
	}
	else {
		return pow(value, a);
	}
}

void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

double cylinder_intersection(double* Ro, double* Rd,
			     double* C, double r) {
  // Step 1. Find the equation for the object you are
  // interested in..  (e.g., cylinder)
  //
  // x^2 + z^2 = r^2
  //
  // Step 2. Parameterize the equation with a center point
  // if needed
  //
  // (x-Cx)^2 + (z-Cz)^2 = r^2
  //
  // Step 3. Substitute the eq for a ray into our object
  // equation.
  //
  // (Rox + t*Rdx - Cx)^2 + (Roz + t*Rdz - Cz)^2 - r^2 = 0
  //
  // Step 4. Solve for t.
  //
  // Step 4a. Rewrite the equation (flatten).
  //
  // -r^2 +
  // t^2 * Rdx^2 +
  // t^2 * Rdz^2 +
  // 2*t * Rox * Rdx -
  // 2*t * Rdx * Cx +
  // 2*t * Roz * Rdz -
  // 2*t * Rdz * Cz +
  // Rox^2 -
  // 2*Rox*Cx +
  // Cx^2 +
  // Roz^2 -
  // 2*Roz*Cz +
  // Cz^2 = 0
  //
  // Steb 4b. Rewrite the equation in terms of t.
  //
  // t^2 * (Rdx^2 + Rdz^2) +
  // t * (2 * (Rox * Rdx - Rdx * Cx + Roz * Rdz - Rdz * Cz)) +
  // Rox^2 - 2*Rox*Cx + Cx^2 + Roz^2 - 2*Roz*Cz + Cz^2 - r^2 = 0
  //
  // Use the quadratic equation to solve for t..
  double a = (sqr(Rd[0]) + sqr(Rd[2]));
  double b = (2 * (Ro[0] * Rd[0] - Rd[0] * C[0] + Ro[2] * Rd[2] - Rd[2] * C[2]));
  double c = sqr(Ro[0]) - 2*Ro[0]*C[0] + sqr(C[0]) + sqr(Ro[2]) - 2*Ro[2]*C[2] + sqr(C[2]) - sqr(r);

  double det = sqr(b) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);
  
  double t0 = (-b - det) / (2*a);
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;
}

double sphere_intersection(double* Ro, double* Rd,
			     double* C, double r) {

	// SAME IDEA AS ABOVE, BUT INCLUDING A Y COMPONENT
	double a = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
	double b = (2*(Ro[0]*Rd[0] - Rd[0]*C[0] + Ro[1]*Rd[1] - Rd[1]*C[1] + Ro[2]*Rd[2] - Rd[2]*C[2]));
	double c = sqr(Ro[0]) - 2*Ro[0]*C[0] + sqr(C[0]) + Ro[1] - 2*Ro[1]*C[1] + sqr(C[1]) + sqr(Ro[2]) - 2*Ro[2]*C[2] + sqr(C[2]) - sqr(r);

  double det = sqr(b) - 4 * a * c;
  if (det < 0) return -1;

  det = sqrt(det);
  
  double t0 = (-b - det) / (2*a);
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;

}

double plane_intersection(double* Ro, double* Rd,
			     double* C, double* N) {

	// USING EQUATIONS FOUND IN THE READING
	double subtract[3];
	subtract[0] = C[0]-Ro[0];
	subtract[1] = C[1]-Ro[1];
	subtract[2] = C[2]-Ro[2];
	double dot1 = N[0]*subtract[0] + N[1]*subtract[1] + N[2]*subtract[2];
	double dot2 = N[0]*Rd[0] + N[1]*Rd[1] + N[2]*Rd[2];
	
	return dot1/dot2;
}

double clamp(double num){
	if (num < 0){
		num = 0;
	}
	else if (num > 1) {
		num = 1;
	}
	return num;
}

int main(int argc, char **argv) {

	// OPEN FILE
	FILE* output;
	output = fopen(argv[4], "wb+");
	// I SHOULD ERROR CHECK HERE FOR THE CORRECT TYPE OF OUTPUT FILE, BUT I'M OUT OF TIME
	// WRITING HEADER INFO
	fprintf(output, "P3\n");
	fprintf(output, "%d %d\n%d\n", atoi(argv[1]), atoi(argv[2]), 1);

	// READING JSON OBJECTS INTO ARRAY  
	read_scene(argv[3]);
	int i = 0;
	double w;
	double h;
	// FINDING CAMERA TO SET WIDTH AND HEIGHT VARIABLES
	while(1){
		if (object_array[i]->kind == 0){
			w = object_array[i]->camera.width;
  			h = object_array[i]->camera.height;
			printf("Camera found and variables set.\n");
			break;
		}
		i++;
	}
  	
	double cx = 0;
  	double cy = 0;

	double d_color[3];
	double s_color[3];

  	int M = atoi(argv[2]);
  	int N = atoi(argv[1]);

  	double pixheight = h / M;
  	double pixwidth = w / N;

	int best;

	// HOLDER VARIABLE FOR CLOSEST OBJECTS COLOR
	double* color;
	// DECREMENTING Y COMPONENT TO FLIP PICTURE 
	for (int y = M; y > 0; y--) {
    	for (int x = 0; x < N; x += 1) {
      		double Ro[3] = {0, 0, 0};
      		// Rd = normalize(P - Ro)
      		double Rd[3] = {
				cx - (w/2) + pixwidth * (x + 0.5),
				cy - (h/2) + pixheight * (y + 0.5),
				1};
      		normalize(Rd);
      		double best_t = INFINITY;
      		for (int i=0; i < obj; i++) {
				//printf("in loop\n");
				double t = 0;
				switch(object_array[i]->kind) {
					//printf("finding best object\n");
					case 0:
	  					// pass, its a camera
	  					break;
					case 1:
						// CHECK INTERSECTION FOR SPHERE
	  					t = sphere_intersection(Ro, Rd,
			    			object_array[i]->sphere.position,
			    			object_array[i]->sphere.radius);
	  						break;
					case 2:
						// CHECK INTERSECTION FOR PLANE
	  					t = plane_intersection(Ro, Rd,
			    			object_array[i]->plane.position,
			    			object_array[i]->plane.normal);
  						break;
  					case 3:
						// pass, its a light
  						break;
					default:
	 					// Horrible error
	 					exit(1);
				}
				t = sqrt(t*t);
				//printf("t = %lf\n", t);
				if (t > 0 && t < best_t) {
					//printf("best t\n");
					best_t = t;
					// SET COLOR TO CLOSEST OBJECTS COLOR
					//color = object_array[i]->color;
					//printf("color stuff\n");
					best = i;
					if (object_array[best]->kind == 1){
						d_color[0] = object_array[best]->sphere.diffuse_color[0];
						d_color[1] = object_array[best]->sphere.diffuse_color[1];
						d_color[2] = object_array[best]->sphere.diffuse_color[2];
					}
					if (object_array[best]->kind == 2){
						d_color[0] = object_array[best]->plane.diffuse_color[0];
						d_color[1] = object_array[best]->plane.diffuse_color[1];
						d_color[2] = object_array[best]->plane.diffuse_color[2];
					}
					if (object_array[best]->kind == 1){
						s_color[0] = object_array[best]->sphere.specular_color[0];
						s_color[1] = object_array[best]->sphere.specular_color[1];
						s_color[2] = object_array[best]->sphere.specular_color[2];
					}
					if (object_array[best]->kind == 2){
						s_color[0] = object_array[best]->plane.specular_color[0];
						s_color[1] = object_array[best]->plane.specular_color[1];
						s_color[2] = object_array[best]->plane.specular_color[2];
					}
					//printf("color stuff end\n");
				}
				//printf("leaving i: %i\n", i);
			}
				
			
			int l = 0;
			double ron[3];
			double rdn[3];
			double distance;
			for (l = 0; l < light; l++){
				//printf("looping lights\n");
				ron[0] = best_t*Rd[0]+Ro[0];
				ron[1] = best_t*Rd[1]+Ro[1];
				ron[2] = best_t*Rd[2]+Ro[2];
				rdn[0] = (lights[l]->light.position[0])-ron[0];
				rdn[1] = (lights[l]->light.position[1])-ron[1];
				rdn[2] = (lights[l]->light.position[2])-ron[2];
				distance = sqrt(rdn[0]*rdn[0] + rdn[1]*rdn[1] + rdn[2]*rdn[2]);
				double closest_shadow = 0;
				int k;
				for (k=0; object_array[k] != 0; k++) {
					//printf("in loop\n");
					double t = 0;
					switch(object_array[k]->kind) {
						case 0:
						// pass, its a camera
		  				break;
						case 1:
						// CHECK INTERSECTION FOR SPHERE
		  				t = sphere_intersection(ron, rdn,
							object_array[k]->sphere.position,
							object_array[k]->sphere.radius);
		  				break;
						case 2:
						// CHECK INTERSECTION FOR PLANE
		  				t = plane_intersection(ron, rdn,
							object_array[k]->plane.position,
							object_array[k]->plane.normal);
		  				break;
						case 3:
						//pass, its another light
		  				break;
						default:
		  				// Horrible error
		  				exit(1);
						}
						if (t > distance || t < 0) {
							continue;
						}
						closest_shadow = 1;
					}
					if (closest_shadow == 0){
						double* N;
		        		double* L;
		        		double* R;
		        		double* V;
		        		double* Kd;
						double* Ks;
					
						if (object_array[best]->kind == 1) {
							N[0] = ron[0]-object_array[best]->sphere.position[0];
							N[1] = ron[1]-object_array[best]->sphere.position[1];
							N[2] = ron[2]-object_array[best]->sphere.position[2];
						}
						else if(object_array[best]->kind == 2) {
							N[0] = object_array[best]->plane.normal[0];
							N[1] = object_array[best]->plane.normal[1];
							N[2] = object_array[best]->plane.normal[2];
						}
						N[0] = sqrt(N[0]*N[0]);
						N[1] = sqrt(N[1]*N[1]);
						N[2] = sqrt(N[2]*N[2]);
						L = rdn;
						R[0] = (2 * dot(L,N) * N[0]) - L[0];
						R[1] = (2 * dot(L,N) * N[1]) - L[1];
						R[2] = (2 * dot(L,N) * N[2]) - L[2];
						V[0] = -Rd[0];
						V[1] = -Rd[1];
						V[2] = -Rd[2];
						Kd = diffuse(d_color, lights[k]->light.color, N, L);
						Ks = specular(s_color, lights[k]->light.color, V, R, N, L);


						double* radial;
						radial[0] = lights[k]->light.radial_a0;
						radial[1] = lights[k]->light.radial_a1;
						radial[2] = lights[k]->light.radial_a2;
						color[0] = frad(radial, distance) * fang(lights[k]->light.angular_a0, lights[k]->light.theta, L, lights[k]->light.direction) * (Kd[0] + Ks[0]);
						color[1] = frad(radial, distance) * fang(lights[k]->light.angular_a0, lights[k]->light.theta, L, lights[k]->light.direction) * (Kd[1] + Ks[1]);
						color[2] = frad(radial, distance) * fang(lights[k]->light.angular_a0, lights[k]->light.theta, L, lights[k]->light.direction) * (Kd[2] + Ks[2]);
					}
				
				}

				// CREATING A PIXEL
				Pixel new;
				//printf("making pixel\n");
    			if (best_t > 0 && best_t != INFINITY) {
					// SETTING PIXELS COLOR TO CLOSEST OBJECTS COLOR
					new.red = color[0];
					new.green = color[1];
					new.blue = color[2];
					// WRITING TO FILE IMMEDIATELY
					printf("plsss");
					fprintf(output, "%i %i %i ", new.red, new.green, new.blue);
    			} else {
					// OTHERWISE ITS AN EMPTY PIXEL, DEFAULT BLACK
					fprintf(output, "0 0 0 ");
    			}
    		}
  		}
  	fclose(output);
  	return 0;
}
