# rev37 Templating Morris Summary

Trajectories r=25; factors=24; template types=4.

## Pendular dose optima (Template_Target_Score vs template_dose)

| Template | Optimal dose (frac of void) | Max score |
|---|---|---|
| capillary_bridge | 0.789 | 0.101 |
| gas | 0.020 | 0.000 |
| rigid | 0.166 | 0.137 |
| surface_weld | 0.240 | 0.239 |

## Top factors for Template_Target_Score by template type

**rigid**: template_dose(1.00,lin), Q_template(0.72,lin), D_template(0.45,lin), C_temp_mass(0.41,lin), C_solid_mass(0.38,lin), Q_colloid(0.27,lin) 

**gas**: C_solid_mass(1.00,lin), Q_template(0.81,cplx), Q_colloid(0.30,lin), v_tip(0.00,cplx), tau(0.00,cplx), P_system(0.00,cplx) 

**surface_weld**: template_dose(1.00,cplx), tau(0.43,cplx), A_molecule(0.27,cplx), C_surfactant(0.25,cplx), MW_surfactant(0.16,cplx), C_plasticizer(0.12,cplx) 

**capillary_bridge**: D_particle(1.00,cplx), template_dose(0.80,cplx), A_molecule(0.75,cplx), MW_surfactant(0.69,cplx), C_surfactant(0.63,cplx), tau(0.18,cplx) 

