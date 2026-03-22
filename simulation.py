from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
from astropy.time import Time
from datetime import datetime
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import rebound
import logging

CSV_FILE_PATH = 'csvs/'
PLOT_FILE_PATH = 'plots/'

G = 4 * np.pi ** 2
M_s = 1.989e30

NAMES = ['sun', 'mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'flyby']
TIME = Time('2026-01-01')

logging.basicConfig(
    level=logging.INFO,
    format="{asctime} | {levelname} | {name} | {lineno} | {funcName} | {message}",
    style="{",
    datefmt='%Y-%m-%d-%H:%M',
)
logger = logging.getLogger(__name__)

# VARY THE LINES BELOW TO ADJUST STELLAR FLYBY PARAMETERS
B = 100
START_DISTANCE = 1000
VELOCITY_AT_INFINITY = 2.1
FLYBY_MASS = 1.0

MASSES = np.array([1.0, 3.301e23/M_s, 4.867e24/M_s, 5.972e24/M_s, 6.417e23/M_s, 1.898e27/M_s, 5.683e26/M_s, 8.618e25/M_s, 1.024e26/M_s, FLYBY_MASS])

def initialise_planets_and_flyby():
    logger.info('Initialising Planets...')
    solar_system_ephemeris.set('de432s')
    positions = np.zeros((10, 3))
    velocities = np.zeros((10, 3))
    for i, body in enumerate(NAMES[:-1]):
        if body == 'sun':
            continue
        pos, vel = get_body_barycentric_posvel(body, TIME)
        pos = pos.xyz.to(u.AU).value
        vel = vel.xyz.to(u.AU / u.year).value
        positions[i] = pos
        velocities[i] = vel

    positions[9] = np.array([START_DISTANCE, B, 0.0])
    velocities[9] = np.array([-VELOCITY_AT_INFINITY, 0.0, 0.0])
    logger.info(f'Planets initialised:\n{positions}\n{velocities}')
    return positions, velocities

def calculate_orbital_elements(r, v):
    mu = 4 * np.pi ** 2

    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)

    specific_orbital_energy = 0.5 * v_mag ** 2 - mu / r_mag

    semi_major_axis = np.inf if specific_orbital_energy >= 0 else -mu / (2 * specific_orbital_energy)

    angular_momentum = np.linalg.norm(np.cross(r, v))

    e_s = 1 + (2 * specific_orbital_energy * angular_momentum ** 2) / mu ** 2
    eccentricity = np.sqrt(max(0, e_s))

    semi_minor_axis = np.inf if eccentricity > 1 else semi_major_axis * np.sqrt(1 - eccentricity ** 2)

    return semi_major_axis, semi_minor_axis, eccentricity, angular_momentum

def initialise_rebound_sim(positions, velocities):
    sim = rebound.Simulation()
    logger.info(f'Initialising simulation: {sim}')
    sim.units = ('yr', 'AU', 'Msun')
    sim.integrator = 'ias15'
    sim.G = 4 * np.pi ** 2

    for i, body in enumerate(NAMES):
        sim.add(m=MASSES[i], x=positions[i, 0], y=positions[i, 1], z=positions[i, 2],
                vx=velocities[i, 0], vy=velocities[i, 1], vz=velocities[i, 2])
        logger.info(f'Added planet {body} with properties mass={MASSES[i]},'
                    f'(x, y, z)=({positions[i, 0]}, {positions[i, 1]}, {positions[i, 2]}),'
                    f'(vx, vy, vz)=({velocities[i, 0]}, {velocities[i, 1]}, {velocities[i, 2]})')

    logger.info(f'Initialised rebound simulation: {sim} with particles {list(sim.particles)}')

    return sim

def run():
    logger.info('Running simulation...')
    fig, ax = plt.subplots(2, 2, constrained_layout=True)
    positions, velocities = initialise_planets_and_flyby()
    sim = initialise_rebound_sim(positions, velocities)

    initial_orbital_elements = {}
    for i, body in enumerate(NAMES):
        if body == 'sun' or body == 'flyby':
            continue
        a, b, e, h = calculate_orbital_elements(positions[i], velocities[i])
        initial_orbital_elements[body] = {'a': a, 'b': b, 'e': e, 'h': h}

    # YOU MAY HAVE TO ADJUST THIS IF YOU ALTER THE START DISTANCE AND VELOCITY
    total_time = 20000

    earth_e = []
    earth_a = []
    earth_b = []
    earth_h = []

    for t in range(total_time):
        if t % 1000 == 0:
            print(f'Step {t} out of {total_time}')
        sim.integrate(t)

        earth = sim.particles[3]
        earth_pos = np.array([earth.x, earth.y, earth.z])
        earth_vel = np.array([earth.vx, earth.vy, earth.vz])

        sun = sim.particles[0]
        sun_pos = np.array([sun.x, sun.y, sun.z])
        sun_vel = np.array([sun.vx, sun.vy, sun.vz])

        rel_pos = earth_pos - sun_pos
        rel_vel = earth_vel - sun_vel

        da, db, de, dh = calculate_orbital_elements(rel_pos, rel_vel)
        earth_e.append(de)
        earth_a.append(da)
        earth_b.append(db)
        earth_h.append(dh)

    final_orbital_elements = {}
    for i, body in enumerate(NAMES):
        if body == 'sun' or body == 'flyby':
            continue
        sun = sim.particles[0]
        curr_body = sim.particles[i]

        pos = np.array([curr_body.x - sun.x, curr_body.y - sun.y, curr_body.z - sun.z])
        vel = np.array([curr_body.vx - sun.vx, curr_body.vy - sun.vy, curr_body.vz - sun.vz])

        a, b, e, h = calculate_orbital_elements(pos, vel)
        final_orbital_elements[body] = {'a': a, 'b': b, 'e': e, 'h': h}

    results = []
    for i, body in enumerate(NAMES):
        if body == 'sun' or body == 'flyby':
            continue
        initial = initial_orbital_elements[body]
        final = final_orbital_elements[body]

        results.append({
            'a1': initial['a'],
            'a2': final['a'],
            'δa': (final['a'] - initial['a']),
            'b1': initial['b'],
            'b2': final['b'],
            'δb': (final['b'] - initial['b']),
            'e1': initial['e'],
            'e2': final['e'],
            'δe': (final['e'] - initial['e']),
            'h1': initial['h'],
            'h2': final['h'],
            'δh': (final['h'] - initial['h']),
        })

    df = pd.DataFrame(results)
    df.to_csv(f'{CSV_FILE_PATH}{datetime.now().strftime('%H_%M')}_{FLYBY_MASS}_{B}.csv', index=False)

    x = np.arange(1, total_time + 1)
    ax[0, 0].plot(x, earth_e); ax[0, 0].set_title('Eccentricity'); ax[0, 0].set_ylabel('e'); ax[0, 0].grid(True)
    ax[0, 1].plot(x, earth_a); ax[0, 1].set_title('Semi-major axis'); ax[0, 1].set_ylabel('a'); ax[0, 1].grid(True)
    ax[1, 0].plot(x, earth_b); ax[1, 0].set_title('Semi-minor axis'); ax[1, 0].set_ylabel('b'); ax[1, 0].grid(True)
    ax[1, 1].plot(x, earth_h); ax[1, 1].set_title('Angular momentum'); ax[1, 1].set_ylabel('h'); ax[1, 1].grid(True)

    plt.savefig(f'{PLOT_FILE_PATH}{datetime.now().strftime('%H_%M')}_{FLYBY_MASS}_{B}.pdf')

    logger.info('Done!')

if __name__ == '__main__':
    run()
