#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


int line = 1;

int width;
int height;

double cam_width;
double cam_height;


typedef struct {
  int kind; // 0 = plane, 1 = sphere, 2 = camera, 3 = light;
  double diffuse_color[3];
  double specular_color[3];
  double reflectivity;
  double refractivity;
  //char sym;
  union {
    struct {
      double position[3]; //center position;
      double normal[3];
    } plane;
    struct {
      double position[3]; //center position;
      double radius;

    } sphere;
    struct {
      double color[3];
      double position[3];
      double direction[3];
      double radial_a2;
      double radial_a1;
      double radial_a0;
      double angular_a1;
      double theta;
    } light;
  };
} Object;


typedef struct shoot_info {
   double t_value;
   int object_index;

} shoot_info;


typedef struct Pixel {

    unsigned char r, g, b;

 }   Pixel;



typedef struct {
    double width;
    double height;

} Camera;

double clamp(double value){
    if(value > 1){
        return 1;
    }else if(value < 0){
        return 0;
    }else{
        return value;
    }
}
typedef double* V3;

static inline void v3_add(V3 a, V3 b, V3 c) {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];


}

static inline void v3_subtract(V3 a, V3 b, V3 c) {

  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

static inline void v3_scale(V3 a, double s, V3 c) {

  c[0] = s * a[0];
  c[1] = s * a[1];
  c[2] = s * a[2];



}

static inline double v3_dot(V3 a, V3 b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline void v3_cross(V3 a, V3 b, V3 c) {

  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];



}

static inline double v3_distance(V3 a, V3 b){
/*
  for(int i = 0; i < 3; i++){
    printf("a value: %f\n", a[i]);
    printf("b value: %f\n", b[i]);
  }
  */
  double dx = b[0] - a[0];
  double dy = b[1] - a[1];
  double dz = b[2] - a[1];
  double c = sqrt(dx*dx + dy*dy + dz*dz);
  //printf("distance value: %d \n", c);
  return c;
}
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
    ;
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Object** read_scene(char* filename, Object** objects) {

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
  int index = 0;
  while (1) {
    //printf("Index: %d \n", index);
    c = fgetc(json);
    if (c == ']') {
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      exit(1);
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

      if (strcmp(value, "camera") == 0) {
            if(index ==0 ){
                index = -1;
            }else{
                index--;
            }
      }else if(strcmp(value, "light") == 0){
          objects[index] = malloc(sizeof(Object));
          objects[index]->light.theta = 0;
          objects[index]->light.radial_a2 = 1;
          objects[index]->light.radial_a1 = 0;
          objects[index]->light.radial_a0 = 0;
          objects[index]->light.color[0] = 1;
          objects[index]->light.color[1] = 1;
          objects[index]->light.color[2] = 1;
          objects[index]->kind = (int) 3;
      } else if (strcmp(value, "sphere") == 0) {
          objects[index] = malloc(sizeof(Object));
          objects[index]->diffuse_color[0] = 0;
          objects[index]->diffuse_color[1] = 0;
          objects[index]->diffuse_color[2] = 0;
          objects[index]->specular_color[0] = 1;
          objects[index]->specular_color[1] = 1;
          objects[index]->specular_color[2] = 1;
          objects[index]->reflectivity = .5;
          objects[index]->refractivity = 1;
          objects[index]->kind = (int) 1;
          //printf("Kind: %d \n", objects[index]->kind);
      } else if (strcmp(value, "plane") == 0) {
          objects[index] = malloc(sizeof(Object));
           objects[index]->diffuse_color[0] = 0;
          objects[index]->diffuse_color[1] = 0;
          objects[index]->diffuse_color[2] = 0;
          objects[index]->specular_color[0] = 1;
          objects[index]->specular_color[1] = 1;
          objects[index]->specular_color[2] = 1;
          objects[index]->reflectivity = .5;
          objects[index]->refractivity = 1;
          objects[index]->kind = (int) 0;
          //printf("Kind: %d \n", objects[index]->kind);
      } else {
	fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
	exit(1);
      }

      skip_ws(json);

    while (1) {
	// , }
	c = next_c(json);
	//printf("C at front of loop: %c\n", c);
	if (c == '}') {
	  // stop parsing this object
	  break;
	} else if (c == ',') {
	  // read another field
	  skip_ws(json);
	  char* key = next_string(json);
	  //printf(key);
	  skip_ws(json);
	  expect_c(json, ':');
	  skip_ws(json);
	  if ((strcmp(key, "width") == 0) ||
	      (strcmp(key, "height") == 0) ||
          (strcmp(key, "radial-a0") == 0) ||
          (strcmp(key, "radial-a1") == 0) ||
          (strcmp(key, "radial-a2") == 0) ||
          (strcmp(key, "angular-a1") == 0) ||
          (strcmp(key, "theta") == 0) ||
          (strcmp(key, "reflectivity") == 0) ||
          (strcmp(key, "refractivity") == 0) ||
	      (strcmp(key, "radius") == 0)) {
	    double value = next_number(json);
	    //printf("VALUE: %f", value);
	    //printf("Value: %f\n", value);
	    if(strcmp(key, "width") == 0){
        cam_width = value;
	    }
	    if(strcmp(key, "height") == 0){
        cam_height = value;
	    }
	    if(strcmp(key, "reflectivity") == 0){
            objects[index]->reflectivity = value;
	    }
	    if(strcmp(key, "refractivity") == 0){
            objects[index]->refractivity = value;
	    }
	    if(strcmp(key, "radius") == 0){
          if(objects[index]->kind == 1){
              objects[index]->sphere.radius = value;
              //printf("radius: %f", value);

          }
	    }else if(strcmp(key, "radial-a0") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.radial_a0 = value;
                //printf("value: %f \n", value);
            }
	    }else if(strcmp(key, "radial-a1") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.radial_a1 = value;
            }
	    }else if(strcmp(key, "radial-a2") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.radial_a2 = value;
            }
	    }else if(strcmp(key, "angular-a1") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.angular_a1 = value;
            }
	    }else if(strcmp(key, "theta") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.theta = value;
            }
	    }
	  }else if ((strcmp(key, "diffuse_color") == 0) ||
             (strcmp(key, "specular_color") == 0) ||
             (strcmp(key, "color") == 0) ||
		     (strcmp(key, "position") == 0) ||
             (strcmp(key, "direction") == 0) ||
		     (strcmp(key, "normal") == 0)) {
	    double* value = next_vector(json);
	    if(strcmp(key, "diffuse_color") == 0){
            objects[index]->diffuse_color[0] = value[0];
            objects[index]->diffuse_color[1] = value[1];
            objects[index]->diffuse_color[2] = value[2];
            //printf("got to color...");
         }else if(strcmp(key, "specular_color") == 0){
            objects[index]->specular_color[0] = value[0];
            objects[index]->specular_color[1] = value[1];
            objects[index]->specular_color[2] = value[2];
         }else if(strcmp(key, "color") == 0){
            if(objects[index]->kind == 3){
                objects[index]->light.color[0] = value[0];
                objects[index]->light.color[1] = value[1];
                objects[index]->light.color[2] = value[2];
            }
         }else if(strcmp(key, "position") == 0){
            if(objects[index]->kind == 0){
                objects[index]->plane.position[0] = value[0];
                objects[index]->plane.position[1] = value[1]; //*-1
                objects[index]->plane.position[2] = value[2];
                //printf("got to position(plane)...");
            }
            else if(objects[index]->kind == 1){
                objects[index]->sphere.position[0] = value[0];
                objects[index]->sphere.position[1] = value[1]; //*-1
                objects[index]->sphere.position[2] = value[2];
            }else if(objects[index]->kind == 3){
                objects[index]->light.position[0] = value[0];
                objects[index]->light.position[1] = value[1];  //*-1
                objects[index]->light.position[2] = value[2];
            }else{
                //this attribute does not belong in this objects
                printf("THIS WENT WRONG");
            }

         }else if(strcmp(key, "normal") == 0){
            if(objects[index]->kind == 0){
                objects[index]->plane.normal[0] = value[0];
                objects[index]->plane.normal[1] = value[1];
                objects[index]->plane.normal[2] = value[2];
            }else
                //this attribute does not belong in this objects
                printf("This went TERRIBLY wrong...");
            }
         else if (strcmp(key, "direction") == 0){

           if(objects[index]->kind == 3){
                objects[index]->light.direction[0] = value[0];
                objects[index]->light.direction[0] = value[0];
                objects[index]->light.direction[0] = value[0];
           }else{
                //this attribute does not belong in this objects
           }
        }
	   }else {
	    fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
		    key, line);
	    //char* value = next_string(json);
	  }

	  skip_ws(json);
	  }else {
	  fprintf(stderr, "Error: Unexpected value on line %d\n", line);
	  exit(1);
	}
      }
      //printf("C is: %c \n", c);
      skip_ws(json);
      //printf("C is: %c \n", c);
      c = next_c(json);
      //printf("C is: %c \n", c);
      if (c == ',') {
        // noop
        skip_ws(json);
      } else if (c == ']') {
        //printf("Closing File....");
        fclose(json);
        objects[index+1] = NULL;
        //printf("Kind: %d \n", objects[index-2]->kind);
        return objects;
      } else {
	fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
	exit(1);
      }
    }
    index++;
  }

}
/////////////////////////////////////////////////////////////////////////////////////


static inline double sqr(double v) {
  return v*v;
}

static inline double to_radians(double degrees){
    return degrees * (M_PI/180.0);
}

static inline void normalize(double* v) {
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

double sphere_intersection(double* Ro, double* Rd,
			     double* C, double r) {

    double a = (sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]));
    double b = (2*(Rd[0]*(Ro[0] - C[0]) + Rd[1]*(Ro[1] - C[1]) + Rd[2]*(Ro[2] - C[2])));
    double c = sqr(Ro[0] - C[0]) + sqr(Ro[1] - C[1]) + sqr(Ro[2] - C[2]) - sqr(r);

    double det = sqr(b) - 4 *a*c;

    //printf("%f ", det);

  if (det < 0){return -1;}

    //printf("test");
    det = sqrt(det);

  double t0 = (-b - det / (2*a));
  if (t0 > 0) return t0;

  double t1 = (-b + det) / (2*a);
  if (t1 > 0) return t1;

  return -1;

}

double plane_intersection(double* Ro, double* Rd, double* L, double* N){
    double D;
    double t;

    for(int i = 0; i <3; i++){
        D += -1*L[i]*N[i];
    }

    t = -(N[0]*Ro[0] + N[1]*Ro[1] + N[2]*Ro[2] + D)/(N[0]*Rd[0] + N[1]*Rd[1] + N[2]*Rd[2]);

    return t;

}

static void reflect_vector(Object* object, double* D, double* R, double* Ron){
    double* N = malloc(sizeof(double)*3);
    if(object->kind == 0){
        N[0] = object->plane.normal[0];
        N[1] = object->plane.normal[1];
        N[2] = object->plane.normal[2];

        v3_scale(N, -1, N);
        normalize(N);
    }else if(object->kind == 1){
        v3_subtract(Ron, object->sphere.position, N); // sphere
        normalize(N);
    }else{
        perror("Object kind not known...");
        exit(1);
    }

    v3_scale(N, 2, R);
    v3_scale(N, v3_dot(R, D), R);
    v3_subtract(R, D, R);
}

shoot_info shoot(double* Ro, double* Rd, Object** scene, int shadow_object){ //set shadow_object to -1 if not testing for shadows
      //printf("running shoot...\n");

      double best_t = INFINITY;
      int current_object = 0;
      int closest_object = -1;
      for (int i=0; scene[i] != 0; i += 1) {
        double t = 0;
        if(scene[i]->kind == 3){ //if light continue
          continue;
        }
    current_object = i;
    if(shadow_object != -1){
        if(shadow_object == current_object){
            //printf("is current object...\n");
            continue;
        }
    }
	switch(scene[i]->kind) {
	case 1:
      //printf("running sphere intersection...\n");
	  t = sphere_intersection(Ro, Rd,
				    scene[i]->sphere.position,
				    scene[i]->sphere.radius);
				    //printf("%f", t);
      //printf("completed sphere intersection...\n");
	  break;
    case 0:
        //printf("running plane intersection...\n");
        t = plane_intersection(Ro, Rd,
                     scene[i]->plane.position,
                     scene[i]->plane.normal);
        //printf("completed plane intersection...\n");
        break;
	default:
	  // Horrible error
	  printf("horrible error");
	  exit(1);
	}


	if (t > 0 && t < best_t){
        best_t = t;
        closest_object = current_object;
	}
      }
      shoot_info result;
      result.t_value = best_t;
      result.object_index = closest_object;
      //final_object = closest_object;

      return result;

}
 double frad(double dl, Object* light){
          double a0 = light->light.radial_a0;
          double a1 = light->light.radial_a1;
          double a2 = light->light.radial_a2;

          //printf("A0 = %f \n", a0);
          //printf("A1 = %f \n", a1);
          //printf("A2 = %f \n", a2);
          double ans =  1/(a0 + a1*dl + a2*sqr(dl));
          //printf("returned %f from frad \n", ans);
          return ans;
      }

      double fang(Object* light, double* Vo){ //Vo = N
          //if not spotlight return 1.0
          double a1 = light->light.angular_a1;
          normalize(light->light.direction);
          if(light->light.theta == 0){
            //printf("returned 1 in fang \n");
            return 1;
          }else if(v3_dot(Vo, light->light.direction) < cos(to_radians(light->light.theta))){ //not in spotlight
            //printf("returned 0 in fang \n");
            return 0;
          }else{
            //printf("returned other value in fang \n");
            return pow(v3_dot(Vo, light->light.direction), a1);
          }

      }

      double diffuse(double light_color, double object_color, double* N, double* L){


            double n_dot_l = v3_dot(N, L);
            double result;
            //printf("N dot L is: %f\n", n_dot_l);
            if(n_dot_l > 0){

                //printf("light color = %f \n", light_color);
                //printf("object folor = %f \n", object_color);
                result =  object_color*light_color*n_dot_l;
                //printf("DIFFUSE VALUE: %f \n", result);
                return result;
            }else{
                return 0;
            }

      }

      double specular(double light_color, double specular_color, double* V, double* R, double* N, double* L){

            double v_dot_r = v3_dot(V, R);
            double n_dot_l = v3_dot(N, L);
            double ns = 20;
            double result;

            if(v_dot_r > 0 && n_dot_l > 0){
                //printf("V DOT R: %f \n", v_dot_r);
                //printf("N DOT L: %f \n", n_dot_l);
                result = specular_color*pow(v_dot_r, ns);
                if(result > 0){
                    //printf("Specular value: %f \n", result);
                }

                return result;
            }else{
                //printf("zero spec value...\n");
                return 0;
                }
      }


void vector_print(double* vector){
   for(int i = 0; i < 3; i++){
      printf("%f\n", vector[i]);
   }
}
double* shade(double* Ro, double* Rd, Object** scene, int current, double best_t, int recursion_level){
    printf("SHADING LEVEL %i: \n", recursion_level);
    double* color = malloc(sizeof(double)*3);
      color[0] = 0;
      color[1] = 0;
      color[2] = 0;

    if(recursion_level == 7){
       return color;
    }else{
     double* Ron = malloc(sizeof(double)*3);
     v3_scale(Rd, best_t, Ron);
     v3_add(Ron, Ro, Ron);

     double* reflected_color = malloc(sizeof(double)*3);
     double* reflection = malloc(sizeof(double)*3);
     reflect_vector(scene[current], Rd, reflection, Ron);
     //printf("IM... GONNA... SHOOT...\n");
     //vector_print(Ron);
     //vector_print(reflection);
     shoot_info new_object = shoot(Ron, reflection, scene, -1);
     //printf("AAAND WE HAVE FIREd!...\n");
     int new_object_index = new_object.object_index;
     double new_t = new_object.t_value;
     if(new_t == INFINITY){
        //printf("infinity and beyond!\n");
        return color;
     }else{
        //printf("GOTTA REFLECT!\n");
        reflected_color = shade(Ron, reflection, scene, new_object_index, new_t, recursion_level +1);
        color[0] += reflected_color[0]*scene[current]->reflectivity;
        color[1] += reflected_color[1]*scene[current]->reflectivity;
        color[2] += reflected_color[2]*scene[current]->reflectivity;
    }



    double* Rdn = malloc(sizeof(double)*3);



      //shadow test
      for(int i = 0; scene[i] != NULL; i++){ //loop through lights
        //if not a light continue...
        if(scene[i]->kind != 3){
            continue;
        }
        //printf("current light a0: %d \n", scene[i]->light.radial_a0);
        //printf("LOOPING THROUGH LIGHTS \n");


        //printf("GOT RON");
        v3_subtract(scene[i]->light.position, Ron, Rdn);
        //printf("GOT RDN");
        normalize(Rdn);

        double distance_to_light = v3_distance(Ron, scene[i]->light.position);

        shoot_info shoot_ray = shoot(Ron, Rdn, scene, current);
        int closest_shadow_object_index = shoot_ray.object_index;
        double new_t = shoot_ray.t_value;
        //printf("shadow object index: %i\n", closest_shadow_object_index);
        Object* closest_shadow_object = scene[closest_shadow_object_index];

        //printf("finished second shoot...\n");

      if(new_t > distance_to_light){
        //printf("no shadow...\n");
        closest_shadow_object = NULL;
      }
      /*
      if(closest_shadow_object_index == current){
        closest_shadow_object = NULL;
      }
      */

         //

      if (closest_shadow_object == NULL) { //not in shadow
        //printf("distance to light: %d\n", distance_to_light);
        //printf("NOT IN SHADOW \n");
        // N, L, R, V
        double* Vobject = malloc(sizeof(double)*3);

        v3_subtract(Ron, scene[i]->light.position, Vobject);
        //printf("v3 subtract done...\n");
        v3_scale(Vobject, -1, Vobject);
        //printf("v3 scale done...\n");
        normalize(Vobject);
        double* N = malloc(sizeof(double)*3);
        double* L = malloc(sizeof(double)*3);
        double* R = malloc(sizeof(double)*3);
        double* V = malloc(sizeof(double)*3);
        printf("%i \n", scene[current]->kind);
        if(scene[current]->kind == 0){
            printf("THIS IS A PLANE \n");
            N[0] = scene[current]->plane.normal[0];
            N[1] = scene[current]->plane.normal[1];
            N[2] = scene[current]->plane.normal[2];
            //v3_scale(N, -1, N);  // plane
            normalize(N);
        }else if(scene[current]->kind == 1){
            printf("THIS IS A SPHERE \n");
            v3_subtract(Ron, scene[current]->sphere.position, N); // sphere
            normalize(N);
        }else{
            printf("error...");
            //something went wrong somewhere :/
        }
        v3_scale(Rdn, -1, L); // light_position - Ron;
        normalize(L);

        v3_scale(N, 2, R);
        v3_scale(N, v3_dot(R, L), R);
        v3_subtract(R, L, R);//reflection of L

        V[0] = Rd[0];
        V[1] = Rd[1];
        V[2] = Rd[2];
        //v3_scale(Rd, -1, V);
        //normalize(V);
        //normalize(R);


         // uses object's diffuse color
         // uses object's specular color


        color[0] += frad(distance_to_light, scene[i]) * fang(scene[i], Vobject) * (diffuse(scene[i]->light.color[0], scene[current]->diffuse_color[0], N, L) + specular(scene[i]->light.color[0], scene[current]->specular_color[0], V, R, N, L));
        color[1] += frad(distance_to_light, scene[i]) * fang(scene[i], Vobject) * (diffuse(scene[i]->light.color[1], scene[current]->diffuse_color[1], N, L) + specular(scene[i]->light.color[1], scene[current]->specular_color[1], V, R, N, L));
        color[2] += frad(distance_to_light, scene[i]) * fang(scene[i], Vobject) * (diffuse(scene[i]->light.color[2], scene[current]->diffuse_color[2], N, L) + specular(scene[i]->light.color[2], scene[current]->specular_color[2], V, R, N, L));



        //free(N);
        //free(L);
        //free(R);
        //free(V);

      }


    }


      }
      //double* reflected_color = shade()

 vector_print(color);
 return color;

}
/*
double* directshade(Object** scene, int objIndex, double* Ron, double* Rd, Object* light){

   double* Rdn = malloc(sizeof(double)*3);

   v3_subtract(light->light.position, Ron, Rdn);
   normalize(Rdn);

   double distance_to_light = v3_distance(Ron, light->light.position);

   double* Vobject = malloc(sizeof(double*3));
   v3_subtract(Ron, light->light.position, Vobject);
   v3_scale(Vobject, -1, Vobject);
   normalize(Vobject);

   double* N = malloc(sizeof(double)*3);
   double* L = malloc(sizeof(double)*3);
   double* R = malloc(sizeof(double)*3);
   double* V = malloc(sizeof(double)*3);

   if(scene[objIndex]->kind == 0){

            N[0] = scene[objIndex]->plane.normal[0];
            N[1] = scene[objIndex]->plane.normal[1];
            N[2] = scene[objIndex]->plane.normal[2];
            v3_scale(N, -1, N);  // plane
            normalize(N);
   }else if(scene[objIndex]->kind == 1){

            v3_subtract(Ron, scene[objIndex]->sphere.position, N); // sphere
            normalize(N);
   }else{
            printf("error...");
            //something went wrong somewhere :/
   }

   v3_scale(Rdn, -1, L);
   normalize(L);

   reflect_vector(scene[objIndex], L, R, Ron);

   V[0] = Rd[0];
   V[1] = Rd[1];
   V[2] = Rd[2];

   color[0] += frad(distance_to_light, light) * fang(light, Vobject) * (diffuse(light->light.color[0], scene[objIndex]->diffuse_color[0], N, L) + specular(light->light.color[0], scene[objIndex]->specular_color[0], V, R, N, L));
   color[1] += frad(distance_to_light, light) * fang(light, Vobject) * (diffuse(light->light.color[1], scene[objIndex]->diffuse_color[1], N, L) + specular(light->light.color[1], scene[objIndex]->specular_color[1], V, R, N, L));
   color[2] += frad(distance_to_light, light) * fang(light, Vobject) * (diffuse(light->light.color[2], scene[objIndex]->diffuse_color[2], N, L) + specular(light->light.color[2], scene[objIndex]->specular_color[2], V, R, N, L));

   return color;

}
*/
Pixel** raytrace(Object** scene){
    Pixel** buffer;


    double cx = 0;
    double cy = 0;
    double h = cam_height;
    double w = cam_width;

    int M = width; //width
    int N = height; //height

    width = M;
    height = N;

    buffer = malloc(sizeof(Pixel*)* M );
    for(int i = 0; i < M; i++){
        buffer[i] = malloc(sizeof(Pixel)*N);
    }
    double pixheight = h / M;
    double pixwidth = w / N;

    int pixcount = 0;

    //FOR EACH PIXEL IN THE SCENE:
    for (int y = 0; y < M; y += 1) {
    for (int x = 0; x < N; x += 1) {
      double Ro[3] = {0, 0, 0};
      double Rd[3] = {
        (cx - (w/2) + pixwidth * (x + 0.5)),
        (cy - (h/2) + pixheight * (y + 0.5)),
        1
      };
      normalize(Rd);




      shoot_info shoot_ray = shoot(Ro, Rd, scene, -1);
      int closest_object_index = shoot_ray.object_index;
      double best_t = shoot_ray.t_value;


       //printf("got t...\n");
       //printf("Object index: %i\n", closest_object_index);
      Object* closest_object = NULL;
      if(closest_object_index > 0){
        closest_object = scene[closest_object_index];
      }

      //set object color
      if(best_t == INFINITY){
        buffer[M-y-1][x].r = 0;
        buffer[M-y-1][x].g = 0;
        buffer[M-y-1][x].b = 0;
        continue;
      }else{ //calculate color of the pixel...
         double* color = malloc(sizeof(double)*3);
         printf("gonna shade the pixel...\n");
         color = shade(Ro, Rd, scene, closest_object_index, best_t, 0);
         buffer[M-y-1][x].r = clamp(color[0])*255;
         buffer[M-y-1][x].g = clamp(color[1])*255;
         buffer[M-y-1][x].b = clamp(color[2])*255;
      }

      }

    }



    return buffer;

}

int writeP6(char* fname, Pixel** buffer){
    FILE* fh;
    fh = fopen(fname, "wb");

    fprintf(fh, "P6\n");
    fprintf(fh, "%d ", width);
    fprintf(fh, "%d\n", height);
    fprintf(fh, "%d\n", 255);


    for(int i = 0; i < width; i++){
        for(int j = 0; j<height; j++){
        unsigned char* rgb;
        rgb = malloc(sizeof(unsigned char)*64);
        rgb[0] = buffer[i][j].r;
        rgb[1] = buffer[i][j].g;
        rgb[2] = buffer[i][j].b;
        fwrite(rgb, 1, 3, fh);
        }


    };

    fclose(fh);


}

int main(int argc, char *argv[]){
    //printf("%d", argc);

    if(argc != 5){
        fprintf(stderr, "Incorrect usage of arguments: width height input.json output.ppm");
        return 1;
    }


    width = atoi(argv[1]);

    height = atoi(argv[2]);

    Object** objects;
    objects = malloc(sizeof(Object*)*144);

    Object** scene = read_scene(argv[3], objects);

    Pixel** buffer = raytrace(scene);

    writeP6(argv[4], buffer);


  return 0;

}
