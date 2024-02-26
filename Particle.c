
/*
Particle type: HALO/CENTRAL/SATELLITE/UNKNOWN




*/
Particle *new_particle(double x, double y, double z, double vx, double vy, double vz, double mass, char *label, int type, int parent_id) {
	Particle *p = malloc(sizeof(*p));
	p->x = x;
	p->y = y;
	p->z = z;
	p->vx = vx;	
	p->vy = vy;
	p->vz = vz;
	p->mass = mass;	
	p->label = label;
	if ((type != HALO) && (type != CENTRAL) && (type != SATELLITE) && (type != UNKNOWN)) {
		printf("particle can be only HALO, CENTRAL, SATELLITE or UNKNOWN\n");
		exit(0);
	} else {
		p->type = type;
	}
	p->parent_id = parent_id;
	return p;
}

int free_particle(Particle *p) {
	free(p->label);
	p->label = NULL;
	free(p);
	p = NULL;
	return 0;
}

int print_particle(Particle *p) {
	char type[15];
	if (p->type == HALO) {
		sprintf(type, "halo");
	} else if (p->type == CENTRAL) {
		sprintf(type, "central");
	} else {
		sprintf(type, "satellite");
	}
	printf("Label: %s, x: %lf, y: %lf, z: %lf, vx: %lf, vy: %lf, vz: %lf, mass: %le, type: %s, parent id: %d\n", p->label, p->x, p->y, p->z, p->vx, p->vy, p->vz, p->mass, type, p->parent_id);
	return 0;
}

Particle *copy_particle(Particle *src) {
	char *label = malloc(20 * sizeof(*label));
	strncpy(label, src->label, 20);
	Particle *p = new_particle(src->x, src->y, src->z, src->vx, src->vy, src->vz, src->mass, label, src->type, src->parent_id);	
	return p;
}

//Ensures PBC given current volume
int PBC(Particle *p) {
	while (p->x < 0.0) {
		p->x += volume_limits[0];		
	}
	while (p->x >= volume_limits[0]) {
		p->x -= volume_limits[0];
	}
	while (p->y < 0.0) {
		p->y += volume_limits[1];
	}
	while (p->y >= volume_limits[1]) {
		p->y -= volume_limits[1];
	}
	while (p->z < 0.0) {
		p->z += volume_limits[2];
	}
	while (p->z >= volume_limits[2]) {
		p->z -= volume_limits[2];
	}
	return 0;
}

int particle_copy_test() {
	Particle *p1 = new_particle(1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, "test particle", HALO, 0);
	Particle *p2 = copy_particle(p1);
	p1->x = 2.0;
	if (p2->x != 1.0) {
		printf("copy failed\n");
		exit(0);
	} else {
		printf("particle test OK \n");
	}
	free(p1);
	free(p2);
	return 0;
}

int particle_tests() {
	particle_copy_test();
	return 0;
}
