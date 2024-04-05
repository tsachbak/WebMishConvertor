using Microsoft.AspNetCore.Mvc;

namespace WebMishConvertor.Controllers
{
    [Route("api/[controller]")]
    public class UploadFileController : Controller
    {
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

            return Ok("File uploaded successfully");
        }
    }
}
