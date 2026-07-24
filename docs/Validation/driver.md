# Validation: driver

Source: `src/morrigan/driver.py`. Tests: `tests/test_driver.py`.

| Test | Pin | Authority |
| --- | --- | --- |
| `test_every_impact_record_is_physically_self_consistent` | Every record's `v_esc` reproduces sqrt(2 G (M_t+M_i)/(R_t+R_i)) rebuilt from the record's own masses and radii; `v_impact` never falls below it | The mutual escape speed (analytical limit) |

The perfect-merger mass sum, the positivity of every extensive quantity, the geometry ranges, the chain handover of masses along a body's impact history, seed determinism, the settings-file error contract, the event-landing timestep, and the closed-form Hill-radius embryo layout are asserted alongside.
