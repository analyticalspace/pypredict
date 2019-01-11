#include <Python.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <termios.h>

#include "predict.h"

struct qth;

// This function was extracted from SingleTrack and shows a number of derived parameters related
// to a satellite observed from a particular point on earth.
//
// NOTE: We don't support uplink, downlink, squint and related parameters.  The related code
//       is convoluted and it's never come up in our usage.  FYI, the 'Edit Transponder Database'
//       menu option is still marked "coming soon" :).  We'll add it back if there's demand.
//
int MakeObservation(double obs_time, struct observation * obs) {
    char geostationary=0, aoshappens=0, decayed=0, visibility=0, sunlit;
    double doppler100=0.0, delay;

    PreCalc(0);
    indx=0;

    if (sat_db.transponders>0)
    {
        PyErr_SetString(PyExc_RuntimeError, "pypredict does not support transponder definition.");
        return -1;
    }

    daynum=obs_time;
    aoshappens=AosHappens(indx);
    geostationary=Geostationary(indx);
    decayed=Decayed(indx,0.0);

    //Calcs
    Calc();
    fk=12756.33*acos(xkmper/(xkmper+sat_alt));

    if (sat_sun_status)
    {
        if (sun_ele<=-12.0 && sat_ele>=0.0) {
            visibility='V';
        } else {
            visibility='D';
        }
    } else {
        visibility='N';
    }
    // gathering power seems much more useful than naked-eye visibility
    sunlit = sat_sun_status;

    doppler100=-100.0e06*((sat_range_rate*1000.0)/299792458.0);
    delay=1000.0*((1000.0*sat_range)/299792458.0);

    //TODO: Seems like FindSun(daynum) should go in here
    FindMoon(daynum);

    //printw(5+tshift,1,"Satellite     Direction     Velocity     Footprint    Altitude     Slant Range");
    //printw(6+tshift,1,"---------     ---------     --------     ---------    --------     -----------");
    //printw(7+tshift,1,"        .            Az           mi            mi          mi              mi");
    //printw(8+tshift,1,"        .            El           km            km          km              km");
    //printw(16+bshift,1,"Eclipse Depth   Orbital Phase   Orbital Model   Squint Angle      AutoTracking");
    //printw(17+bshift,1,"-------------   -------------   -------------   ------------      ------------");

    obs->norad_id = sat.catnum;
    strncpy(&(obs->name), &(sat.name), sizeof(obs->name));
    obs->epoch = (daynum+3651.0)*(86400.0); //See daynum=((start/86400.0)-3651.0);
    obs->latitude = sat_lat;
    obs->longitude = sat_lon;
    obs->azimuth = sat_azi;
    obs->elevation = sat_ele;
    obs->orbital_velocity = 3600.0*sat_vel;
    obs->footprint = fk;
    obs->altitude = sat_alt;
    obs->slant_range = sat_range;
    obs->eclipse_depth = eclipse_depth/deg2rad;
    obs->orbital_phase = 256.0*(phase/twopi);
    strncpy(&(obs->orbital_model), &(ephem), sizeof(obs->orbital_model));
    obs->visibility = visibility;
    obs->sunlit = sunlit;
    obs->orbit = rv;
    obs->geostationary = geostationary;
    obs->has_aos = aoshappens;
    obs->decayed = decayed;
    obs->doppler = doppler100;
    obs->inclination = sat.incl;
    obs->raan = sat.raan;
    obs->argument_perigee = sat.argper;
    obs->perigee = sat.perigee;
    obs->apogee = sat.apogee;
    obs->eccentricity = sat.eccn;
    obs->semi_major_axis = sat.aodp;
    obs->mean_anomaly = sat.meanan;
    obs->bstar = sat.bstar;
    obs->x = sat.x;
    obs->y = sat.y;
    obs->z = sat.z;
    obs->x_vel = sat.x_vel;
    obs->y_vel = sat.y_vel;
    obs->z_vel = sat.z_vel;
    return 0;
}

void PrintObservation(struct observation * obs) {
    printf("NORAD_ID           %ld\n", obs->norad_id);
    printf("Name               %s\n", obs->name);
    printf("Date(epoch)        %f\n", obs->epoch);
    printf("Latitude(N)        %f\n", obs->latitude);
    printf("Longitude(E)       %f\n", obs->longitude);
    printf("Azimuth            %f\n", obs->azimuth);
    printf("Elevation          %f\n", obs->elevation);
    printf("Velocity(km/h)     %f\n", obs->orbital_velocity);
    printf("Footprint          %f\n", obs->footprint);
    printf("Altitude(km)       %f\n", obs->altitude);
    printf("Slant_Range(km)    %f\n", obs->slant_range);
    printf("Eclipse_Depth      %f\n", obs->eclipse_depth);
    printf("Orbital_Phase      %f\n", obs->orbital_phase);
    printf("Orbital_Model      %s\n", obs->orbital_model);
    printf("Visibility         %c\n", obs->visibility);
    printf("Sunlight           %d\n", obs->sunlit);
    printf("Orbit_Number       %ld\n",obs->orbit);
    printf("Geostationary      %d\n", obs->geostationary);
    printf("AOS_Happens        %d\n", obs->has_aos);
    printf("Decayed            %d\n", obs->decayed);
    printf("Doppler            %f\n", obs->doppler);
    printf("Inclination        %f\n", obs->inclination);
    printf("RAAN               %f\n", obs->raan);
    printf("Arg Perigee        %f\n", obs->argument_perigee);
    printf("Perigee            %f\n", obs->perigee);
    printf("Apogee             %f\n", obs->apogee);
    printf("Eccentricity       %f\n", obs->eccentricity);
    printf("Semi Major Axis(km)%f\n", obs->semi_major_axis);
    printf("Mean Anomaly       %f\n", obs->mean_anomaly);
    printf("BStar              %f\n", obs->bstar);
    printf("X                  %f\n", obs->x);
    printf("Y                  %f\n", obs->y);
    printf("Z                  %f\n", obs->z);
    printf("X Velocity         %f\n", obs->x_vel);
    printf("Y Velocity         %f\n", obs->y_vel);
    printf("Z Velocity         %f\n", obs->z_vel);
}

PyObject * PythonifyObservation(observation * obs) {
    //TODO: Add reference count?
#if PY_MAJOR_VERSION >= 3
    return Py_BuildValue("{s:l,s:s,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:s,s:C,s:i,s:l,s:i,s:i,s:i,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
#else
    return Py_BuildValue("{s:l,s:s,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:s,s:c,s:i,s:l,s:i,s:i,s:i,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
#endif
        "norad_id", obs->norad_id,
        "name", obs->name,
        "epoch", obs->epoch,
        "latitude", obs->latitude,
        "longitude", obs->longitude,
        "azimuth", obs->azimuth,
        "elevation", obs->elevation,
        "orbital_velocity", obs->orbital_velocity,
        "footprint", obs->footprint,
        "altitude", obs->altitude,
        "slant_range", obs->slant_range,
        "eclipse_depth", obs->eclipse_depth,
        "orbital_phase", obs->orbital_phase,
        "orbital_model", obs->orbital_model,
        "visibility", obs->visibility,
        "sunlit", obs->sunlit,
        "orbit", obs->orbit,
        "geostationary", obs->geostationary,
        "has_aos", obs->has_aos,
        "decayed", obs->decayed,
        "doppler", obs->doppler,
        "inclination", obs->inclination,
        "raan", obs->raan,
        "arg_perigee", obs->argument_perigee,
        "perigee", obs->perigee,
        "apogee", obs->apogee,
        "eccentricity", obs->eccentricity,
        "semi_major_axis", obs->semi_major_axis,
        "mean_anomaly", obs->mean_anomaly,
        "bstar", obs->bstar,
        "x", obs->x,
        "y", obs->y,
        "z", obs->z,
        "x_vel", obs->x_vel,
        "y_vel", obs->y_vel,
        "z_vel", obs->z_vel
    );
}
char load(PyObject *args) {
    //TODO: Not threadsafe, detect and raise warning?
    int x;
    char *env=NULL;

    /* Set up translation table for computing TLE checksums */
    for (x=0; x<=255; val[x]=0, x++);
    for (x='0'; x<='9'; val[x]=x-'0', x++);

    val['-']=1;

    double epoch;
    const char *tle0, *tle1, *tle2;

    if (!PyArg_ParseTuple(args, "(sss)|d(ddi)",
        &tle0, &tle1, &tle2, &epoch, &qth.stnlat, &qth.stnlong, &qth.stnalt))
    {
        // PyArg_ParseTuple will set appropriate exception string
        return -1;
    };

    if (ReadTLE(tle0,tle1,tle2) != 0)
    {
        PyErr_SetString(PyExc_RuntimeError, "Unable to process TLE");
        return -1;
    }

    // If time isn't set, use current time.
    if (PyObject_Length(args) < 2)
    {
        daynum=CurrentDaynum();
    }
    else
    {
        daynum=((epoch/86400.0)-3651.0);
    }

    // If we haven't already set groundstation location, use predict's default.
    if (PyObject_Length(args) < 3)
    {
        FILE *fd;
        env=getenv("HOME");
        sprintf(qthfile,"%s/.predict/predict.qth",env);
        fd=fopen(qthfile,"r");
        if (fd!=NULL)
        {
            fgets(qth.callsign,16,fd);
            qth.callsign[strlen(qth.callsign)-1]=0;
            fscanf(fd,"%lf", &qth.stnlat);
            fscanf(fd,"%lf", &qth.stnlong);
            fscanf(fd,"%d", &qth.stnalt);
            fclose(fd);
        } else {
            PyErr_SetString(PyExc_RuntimeError, "QTH file could not be loaded.");
            return -1;
        }
    }
    obs_geodetic.lat=qth.stnlat*deg2rad;
    obs_geodetic.lon=-qth.stnlong*deg2rad;
    obs_geodetic.alt=((double)qth.stnalt)/1000.0;
    obs_geodetic.theta=0.0;

    return 0;
}

static PyObject* quick_find(PyObject* self, PyObject *args)
{
    struct observation obs = { 0 };

    if (load(args) != 0 || MakeObservation(daynum, &obs) != 0)
    {
        // load or MakeObservation will set appropriate exceptions if either fails.
        return NULL;
    }

    return PythonifyObservation(&obs);
}

static char quick_find_docs[] =
    "quick_find((tle_line0, tle_line1, tle_line2), time, (gs_lat, gs_lon, gs_alt))\n";

static PyObject* quick_predict(PyObject* self, PyObject *args)
{
    double now;
    int lastel=0;
    char errbuff[100];
    observation obs = { 0 };

    PyObject* transit = PyList_New(0);

    if (transit == NULL)
    {
        goto cleanup_and_raise_exception;
    }

    now=CurrentDaynum();

    if (load(args) != 0)
    {
        // load will set the appropriate exception string if it fails.
        goto cleanup_and_raise_exception;
    }

    //TODO: Seems like this should be based on the freshness of the TLE, not wall clock.
    if ((daynum<now-365.0) || (daynum>now+365.0))
    {
        sprintf(errbuff, "time %s too far from present\n", Daynum2String(daynum));
        PyErr_SetString(PyExc_RuntimeError, errbuff);
        goto cleanup_and_raise_exception;
    }

    PreCalc(0);
    Calc();

    if (MakeObservation(daynum, &obs) != 0)
    {
        // MakeObservation will set appropriate exception string
        goto cleanup_and_raise_exception;
    }

    if (!AosHappens(0))
    {
        sprintf(errbuff, "%lu does not rise above horizon. No AOS.\n", sat.catnum);
        PyErr_SetString(PyExc_RuntimeError, errbuff);
        goto cleanup_and_raise_exception;
    }

    if (Geostationary(0)!=0)
    {
        sprintf(errbuff, "%lu is geostationary.  Does not transit.\n", sat.catnum);
        PyErr_SetString(PyExc_RuntimeError, errbuff);
        goto cleanup_and_raise_exception;
    }

    if (Decayed(indx,daynum)!=0)
    {
        sprintf(errbuff, "%lu has decayed. Cannot calculate transit.\n", sat.catnum);
        PyErr_SetString(PyExc_RuntimeError, errbuff);
        goto cleanup_and_raise_exception;
    }

    /* Make Predictions */
    daynum=FindAOS();

    /* Construct the pass */
    PyObject * py_obs;

    while (iel>=0)
    {
        if (MakeObservation(daynum, &obs) != 0)
        {
            //MakeObservation will set appropriate exception string
            goto cleanup_and_raise_exception;
        }

        py_obs = PythonifyObservation(&obs);

        if (py_obs == NULL) {
            sprintf(errbuff, "Failed to PythonifyObservation. Cannot calculate transit.\n");
            PyErr_SetString(PyExc_RuntimeError, errbuff);
            goto cleanup_and_raise_exception;
        }

        if (PyList_Append(transit, py_obs) != 0)
        {
            goto cleanup_and_raise_exception;
        }

        lastel=iel;
        daynum+=cos((sat_ele-1.0)*deg2rad)*sqrt(sat_alt)/25000.0;
        Calc();
    }

    if (lastel!=0)
    {
        daynum=FindLOS();

        //TODO: FindLOS can fail.  Detect and log warning that transit end is approximate.
        if (daynum > 0) {
            Calc();

            if (MakeObservation(daynum, &obs) != 0)
            {
                goto cleanup_and_raise_exception;
            }

            py_obs = PythonifyObservation(&obs);

            if (py_obs == NULL) {
                sprintf(errbuff, "Failed to PythonifyObservation. Cannot calculate transit.\n");
                goto cleanup_and_raise_exception;
            }

            if (PyList_Append(transit, py_obs) != 0)
            {
                goto cleanup_and_raise_exception;
            }
        }
    }

    return transit;

cleanup_and_raise_exception:
    Py_XDECREF(transit);
    return NULL;
}

static char quick_predict_docs[] =
    "quick_predict((tle_line0, tle_line1, tle_line2), time, (gs_lat, gs_lon, gs_alt))\n";

static PyMethodDef pypredict_funcs[] = {
    {"quick_find"   , (PyCFunction)quick_find   , METH_VARARGS, quick_find_docs},
    {"quick_predict", (PyCFunction)quick_predict, METH_VARARGS, quick_predict_docs},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
    static int cpredict_traverse(PyObject *m, visitproc visit, void *arg) {
        return 0;
    }

    static int cpredict_clear(PyObject *m) {
        return 0;
    }

    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cpredict",
        "Python port of the predict open source satellite tracking library",
        0,
        pypredict_funcs,
        NULL,
        cpredict_traverse,
        cpredict_clear,
        NULL
    };

    PyMODINIT_FUNC
    PyInit_cpredict(void)
    {
        PyObject *module = PyModule_Create(&moduledef);

        if (module == NULL) {
            fprintf(stderr, "ERROR: Unable to initialize python module 'cpredict'\n");
            return NULL;
        }

        return module;
    }
#else
    void initcpredict(void)
    {
        PyObject *m;
        m = Py_InitModule3("cpredict", pypredict_funcs,
                        "Python port of the predict open source satellite tracking library");

        if (m == NULL) {
            fprintf(stderr, "ERROR: Unable to initialize python module 'cpredict'\n");
        }
    }
#endif

