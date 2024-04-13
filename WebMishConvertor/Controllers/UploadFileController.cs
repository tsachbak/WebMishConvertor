using Microsoft.AspNetCore.Mvc;
using OfficeOpenXml;
using System.Collections.ObjectModel;
using System.ComponentModel;
using WebMishConvertor.Models;

namespace WebMishConvertor.Controllers
{
    [Route("api/[controller]")]
    public class UploadFileController : Controller
    {
        private List<double> m_latList;
        private List<double> m_lonList;
        private Dictionary<int, SingleDotPosition> m_dots;

        // GET: Displays the upload form
        [HttpGet]
        public IActionResult Index()
        {
            return View();
        }

        //POST: Handles the file upload
        [HttpPost("uploadfile")]
        public async Task<IActionResult> UploadFile(IFormFile file)
        {
            if (file == null || file.Length == 0) 
            {
                return BadRequest("No file uploaded or file is empty");
            }

            if (await TakeDataFromFile(file))
            {
                SaveDataToDictionary();
                return Ok(Dots);
            }

            return BadRequest("Error processing file");
        }

        private async Task<bool> TakeDataFromFile(IFormFile file)
        {
            m_latList = new List<double>();
            m_lonList = new List<double>();

            ExcelPackage.LicenseContext = OfficeOpenXml.LicenseContext.NonCommercial;

            try
            {
                using (var stram = new MemoryStream())
                {
                    await file.CopyToAsync(stram);
                    stram.Position = 0;

                    using (var package = new ExcelPackage(stram))
                    {
                        ExcelWorksheet ws = package.Workbook.Worksheets[0];

                        for (int row = 1; row <= ws.Dimension.Rows; ++row)
                        {
                            for (int col = 1; col <= ws.Dimension.Columns; ++col)
                            {
                                object cellVal = ws.Cells[row, col].Value;
                                if (col == 1)
                                {
                                    if (cellVal != null && double.TryParse(cellVal.ToString(), out double dblVal))
                                    {
                                        m_latList.Add(dblVal);
                                    }
                                }
                                else if (col == 2)
                                {
                                    if (cellVal != null && double.TryParse(cellVal.ToString(), out double dblVal))
                                    {
                                        m_lonList.Add(dblVal);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                return false;
            }

            return true;
        }

        private void SaveDataToDictionary()
        {
            m_dots = new Dictionary<int, SingleDotPosition>();

            for (int i = 0; i < m_latList.Count; i++)
            {
                m_dots.Add(i + 1, new SingleDotPosition(m_latList[i], m_lonList[i]));
            }
        }

        private Dictionary<int, SingleDotPosition> Dots
        {
            get { return m_dots; }
        }
    }
}
