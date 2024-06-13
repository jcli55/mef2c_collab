import loompy

loom_files = [
"/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/M2C_21d_f_con/M2C_21d_f_con.loom",
"/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/M2C_21d_f_ko/M2C_21d_f_ko.loom",
"/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/M2C_21d_m_con/M2C_21d_m_con.loom",
"/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/M2C_21d_m_ko_repeat_final/M2C_21d_m_ko_repeat_final.loom",
"/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/M2C_29d_m_con_final/M2C_29d_m_con_final.loom",
"/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/M2C_29d_m_ko_final/M2C_29d_m_ko_final.loom"
]
loompy.combine(
    loom_files,
    "/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/merged.loom",
)
